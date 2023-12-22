#include <Eigen/Core>
#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <vector>
#include <thread>

#include "Bank.h"
#include "ComplexFunction.h"
#include "omp_utils.h"
#include "parallel.h"
#include "trees/FunctionTree.h"

#ifdef MRCPP_HAS_OMP
#define mrcpp_get_max_threads() omp_get_max_threads()
#define mrcpp_get_num_procs() omp_get_num_procs()/2
#define mrcpp_set_dynamic(n) omp_set_dynamic(n)
#else
#define mrcpp_get_max_threads() 1
#define mrcpp_get_num_procs() 1
#define mrcpp_get_num_threads() 1
#define mrcpp_get_thread_num() 0
#define mrcpp_set_dynamic(n)
#endif

using mrcpp::Printer;
using namespace std;

namespace mrcpp {

namespace omp {

int n_threads = mrcpp_get_max_threads();

} // namespace omp

using namespace Eigen;

Bank dataBank;

namespace mpi {

bool numerically_exact = false;
int shared_memory_size = 1000;

// these parameters set by initialize()
int world_size = 1;
int world_rank = 0;
int wrk_size = 1;
int wrk_rank = 0;
int share_size = 1;
int share_rank = 0;
int sh_group_rank = 0;
int is_bank = 0;
int is_centralbank = 0;
int is_bankclient = 1;
int is_bankmaster = 0; // only one bankmaster is_bankmaster
int bank_size = 0;
int omp_threads = -1; // can be set to force number of threads
int tot_bank_size = 0; // size of bank, including the task manager
int max_tag = 0;       // max value allowed by MPI
vector<int> bankmaster;
int task_bank = -1; // world rank of the task manager

MPI_Comm comm_wrk;
MPI_Comm comm_share;
MPI_Comm comm_sh_group;
MPI_Comm comm_bank;

} // namespace mpi

int id_shift; // to ensure that nodes, orbitals and functions do not collide

extern int metadata_block[3]; // can add more metadata in future
extern int const size_metadata = 3;

void mpi::initialize() {
    Eigen::setNbThreads(1);
    mrcpp_set_dynamic(0);

#ifdef MRCPP_HAS_MPI
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::world_rank);

    // divide the world into groups
    // each group has its own group communicator definition

    // define independent group of MPI processes, that are not part of comm_wrk
    // for now the new group does not include comm_share
    mpi::comm_bank = MPI_COMM_WORLD; // clients and master
    MPI_Comm comm_remainder;         // clients only

    // set bank_size automatically if not defined by user
    if (mpi::world_size < 2) {
        mpi::bank_size = 0;
    } else if (mpi::bank_size < 0) {
        mpi::bank_size = max(mpi::world_size / 3, 1);
    }
    if (mpi::world_size - mpi::bank_size < 1) MSG_ABORT("No MPI ranks left for working!");
    if (mpi::bank_size < 1 and mpi::world_size > 1) MSG_ABORT("Bank size must be at least one when using MPI!");

    mpi::bankmaster.resize(mpi::bank_size);
    for (int i = 0; i < mpi::bank_size; i++) {
        mpi::bankmaster[i] = mpi::world_size - i - 1; // rank of the bankmasters
    }
    if (mpi::world_rank < mpi::world_size - mpi::bank_size) {
        // everything which is left
        mpi::is_bank = 0;
        mpi::is_centralbank = 0;
        mpi::is_bankclient = 1;
    } else {
        // special group of centralbankmasters
        mpi::is_bank = 1;
        mpi::is_centralbank = 1;
        mpi::is_bankclient = 0;
        if (mpi::world_rank == mpi::world_size - mpi::bank_size) mpi::is_bankmaster = 1;
    }
    MPI_Comm_split(MPI_COMM_WORLD, mpi::is_bankclient, mpi::world_rank, &comm_remainder);

    // split world into groups that can share memory
    MPI_Comm_split_type(comm_remainder, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &mpi::comm_share);

    MPI_Comm_rank(mpi::comm_share, &mpi::share_rank);
    MPI_Comm_size(mpi::comm_share, &mpi::share_size);

    // define a rank of the group
    MPI_Comm_split(comm_remainder, mpi::share_rank, mpi::world_rank, &mpi::comm_sh_group);
    // mpiShRank is color (same color->in same group)
    // MPI_worldrank is key (orders rank within the groups)

    // we define a new orbital rank, so that the orbitals within
    // a shared memory group, have consecutive ranks
    MPI_Comm_rank(mpi::comm_sh_group, &mpi::sh_group_rank);

    mpi::wrk_rank = mpi::share_rank + mpi::sh_group_rank * mpi::world_size;
    MPI_Comm_split(comm_remainder, 0, mpi::wrk_rank, &mpi::comm_wrk);
    // 0 is color (same color->in same group)
    // mpiOrbRank is key (orders rank in the group)

    MPI_Comm_rank(mpi::comm_wrk, &mpi::wrk_rank);
    MPI_Comm_size(mpi::comm_wrk, &mpi::wrk_size);

    // if bank_size is large enough, we reserve one as "task manager"
    mpi::tot_bank_size = mpi::bank_size;
    if (mpi::bank_size <= 2 and mpi::bank_size > 0) {
        // use the first bank as task manager
        mpi::task_bank = mpi::bankmaster[0];
    } else if (mpi::bank_size > 1) {
        // reserve one bank for task management only
        mpi::bank_size--;
        mpi::task_bank = mpi::bankmaster[mpi::bank_size]; // the last rank is reserved as task manager
    }

    // determine the maximum value alowed for mpi tags
    void *val;
    int flag;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &val, &flag); // max value allowed by MPI for tags
    max_tag = *(int *)val / 2;
    id_shift = max_tag / 2; // half is reserved for non orbital.

    // determine the number of threads we can assign to each mpi worker.
    // mrcpp_get_num_procs is total number of hardware logical threads accessible by this mpi
    // We assume that half of them are physical cores.
    // mrcpp_get_max_threads is OMP_NUM_THREADS (environment variable).
    // omp_threads_available is the total number of logical threads available on this compute-node
    // We assume that half of them are physical cores.
    //
    // six conditions should be satisfied:
    // 1) no one use more than mrcpp_get_num_procs()/2
    // 2) NOT ENFORCED: no one use more than mrcpp_get_max_threads, as defined by rank 0
    // 3) the total number of threads used on the compute-node must not exceed omp_threads_available/2
    // 4) Bank needs only one thread
    // 5) workers need as many threads as possible
    // 6) at least one thread

    MPI_Comm comm_share_world;//all that share the memory
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_share_world);

    int n_bank_thisnode; //number of banks on this node
    MPI_Allreduce(&is_bank, &n_bank_thisnode, 1, MPI_INT, MPI_SUM, comm_share_world);
    int n_wrk_thisnode; //number of workers on this node
    MPI_Allreduce(&is_bankclient, &n_wrk_thisnode, 1, MPI_INT, MPI_SUM, comm_share_world);

    int omp_threads_available = thread::hardware_concurrency();
    int nthreads = 1;
    if (is_bankclient) nthreads = (omp_threads_available/2-n_bank_thisnode)/n_wrk_thisnode; // 3) and 5)

    // do not exceed total number of cores accessible (assumed to be half the number of logical threads)
    nthreads = min(nthreads, mrcpp_get_num_procs()); // 1)

    // NB: we do not use OMP_NUM_THREADS. Use all cores accessible. Could change this in the future
    // if OMP_NUM_THREADS is set, do not exceed
    // we enforce that all compute nodes use the same OMP_NUM_THREADS. Rank 0 decides.
    /* int my_OMP_NUM_THREADS = mrcpp_get_max_threads();
    MPI_Bcast(&my_OMP_NUM_THREADS, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (my_OMP_NUM_THREADS > 0) nthreads = min(nthreads, my_OMP_NUM_THREADS); // 2)
    */

    nthreads = max(1, nthreads); // 6)

    if (is_bank) nthreads = 1; // 4)

    //cout<<world_rank<<" found "<<omp_threads_available<<" available threads. omp:"<<omp_get_num_procs()<<" "<<omp_get_max_threads()<<" "<<mrcpp::omp::n_threads<<" On this node: "<<n_bank_thisnode<<" banks "<<n_wrk_thisnode<<" workers"<<" "<<nthreads<<" is bank "<<is_bank<<endl;

    if (omp_threads > 0) {
        if (omp_threads != nthreads and world_rank == 0) {
            cout<<"Warning: recommended number of threads is "<<nthreads<<endl;
            cout<<"setting number of threads to omp_threads, "<<omp_threads<<endl;
        }
        nthreads = omp_threads;
    }

    omp::n_threads = nthreads;
    mrcpp::set_max_threads(nthreads);

    if (mpi::is_bank) {
        // bank is open until end of program
        if (mpi::is_centralbank) { dataBank.open(); }
        mpi::finalize();
        exit(EXIT_SUCCESS);
    }
#else
    mpi::bank_size = 0;
    mrcpp::set_max_threads(omp::n_threads);
#endif
}

void mpi::finalize() {
#ifdef MRCPP_HAS_MPI
    if (mpi::bank_size > 0 and mpi::grand_master()) {
        println(4, " max data in bank " << dataBank.get_maxtotalsize() << " MB ");
        dataBank.close();
    }
    MPI_Barrier(MPI_COMM_WORLD); // to ensure everybody got here
    MPI_Finalize();
#endif
}

void mpi::barrier(MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    MPI_Barrier(comm);
#endif
}

/*********************************
 * Orbital related MPI functions *
 *********************************/

bool mpi::grand_master() {
    return (mpi::world_rank == 0 and is_bankclient) ? true : false;
}

bool mpi::share_master() {
    return (mpi::share_rank == 0) ? true : false;
}

/** @brief Test if orbital belongs to this MPI rank (or is common)*/
bool mpi::my_orb(int j) {
    return ((j) % mpi::wrk_size == mpi::wrk_rank) ? true : false;
}

/** @brief Test if orbital belongs to this MPI rank (or is common)*/
bool mpi::my_orb(ComplexFunction orbj) {
    return my_orb(orbj.getRank());
}

/** @brief Free all function pointers not belonging to this MPI rank */
void mpi::free_foreign(MPI_FuncVector &Phi) {
    for (ComplexFunction &i : Phi) {
        if (not mpi::my_orb(i)) i.free(NUMBER::Total);
    }
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(IntVector &vec, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(DoubleVector &vec, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(ComplexVector &vec, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(IntMatrix &mat, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

// send a function with MPI
void mpi::send_function(ComplexFunction &func, int dst, int tag, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    if (func.isShared()) MSG_WARN("Sending a shared function is not recommended");
    FunctionData &funcinfo = func.getFunctionData();
    MPI_Send(&funcinfo, sizeof(FunctionData), MPI_BYTE, dst, 0, comm);
    if (func.hasReal()) mrcpp::send_tree(func.real(), dst, tag, comm, funcinfo.real_size);
    if (func.hasImag()) mrcpp::send_tree(func.imag(), dst, tag + 10000, comm, funcinfo.imag_size);
#endif
}

// receive a function with MPI
void mpi::recv_function(ComplexFunction &func, int src, int tag, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    if (func.isShared()) MSG_WARN("Receiving a shared function is not recommended");
    MPI_Status status;

    FunctionData &funcinfo = func.getFunctionData();
    MPI_Recv(&funcinfo, sizeof(FunctionData), MPI_BYTE, src, 0, comm, &status);
    if (funcinfo.real_size > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (not func.hasReal()) func.alloc(NUMBER::Real);
        mrcpp::recv_tree(func.real(), src, tag, comm, funcinfo.real_size);
    }

    if (funcinfo.imag_size > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (not func.hasImag()) func.alloc(NUMBER::Imag);
        mrcpp::recv_tree(func.imag(), src, tag + 10000, comm, funcinfo.imag_size);
    }
#endif
}

/** Update a shared function after it has been changed by one of the MPI ranks. */
void mpi::share_function(ComplexFunction &func, int src, int tag, MPI_Comm comm) {
    if (func.isShared()) {
#ifdef MRCPP_HAS_MPI
        if (func.hasReal()) mrcpp::share_tree(func.real(), src, tag, comm);
        if (func.hasImag()) mrcpp::share_tree(func.imag(), src, 2 * tag, comm);
#endif
    }
}

/** @brief Add all mpi function into rank zero */
void mpi::reduce_function(double prec, ComplexFunction &func, MPI_Comm comm) {
/* 1) Each odd rank send to the left rank
   2) All odd ranks are "deleted" (can exit routine)
   3) new "effective" ranks are defined within the non-deleted ranks
      effective rank = rank/fac , where fac are powers of 2
   4) repeat
 */
#ifdef MRCPP_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) {
        if ((comm_rank / fac) % 2 == 0) {
            // receive
            int src = comm_rank + fac;
            if (src < comm_size) {
                ComplexFunction func_i(false);
                int tag = 3333 + src;
                mpi::recv_function(func_i, src, tag, comm);
                func.add(1.0, func_i); // add in place using union grid
                func.crop(prec);
            }
        }
        if ((comm_rank / fac) % 2 == 1) {
            // send
            int dest = comm_rank - fac;
            if (dest >= 0) {
                int tag = 3333 + comm_rank;
                mpi::send_function(func, dest, tag, comm);
                break; // once data is sent we are done
            }
        }
        fac *= 2;
    }
    MPI_Barrier(comm);
#endif
}

/** @brief make union tree and send into rank zero */
void mpi::reduce_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, MPI_Comm comm) {
/* 1) Each odd rank send to the left rank
   2) All odd ranks are "deleted" (can exit routine)
   3) new "effective" ranks are defined within the non-deleted ranks
      effective rank = rank/fac , where fac are powers of 2
   4) repeat
 */
#ifdef MRCPP_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) {
        if ((comm_rank / fac) % 2 == 0) {
            // receive
            int src = comm_rank + fac;
            if (src < comm_size) {
                int tag = 3333 + src;
                mrcpp::FunctionTree<3> tree_i(tree.getMRA());
                mrcpp::recv_tree(tree_i, src, tag, comm, -1, false);
                tree.appendTreeNoCoeff(tree_i); // make union grid
            }
        }
        if ((comm_rank / fac) % 2 == 1) {
            // send
            int dest = comm_rank - fac;
            if (dest >= 0) {
                int tag = 3333 + comm_rank;
                mrcpp::send_tree(tree, dest, tag, comm, -1, false);
                break; // once data is sent we are done
            }
        }
        fac *= 2;
    }
    MPI_Barrier(comm);
#endif
}

/** @brief make union tree without coeff and send to all
 *  Include both real and imaginary parts
 */
void mpi::allreduce_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, vector<ComplexFunction> &Phi, MPI_Comm comm) {
    /* 1) make union grid of own orbitals
       2) make union grid with others orbitals (sent to rank zero)
       3) rank zero broadcast func to everybody
     */

    int N = Phi.size();
    for (int j = 0; j < N; j++) {
        if (not mpi::my_orb(j)) continue;
        if (Phi[j].hasReal()) tree.appendTreeNoCoeff(Phi[j].real());
        if (Phi[j].hasImag()) tree.appendTreeNoCoeff(Phi[j].imag());
    }
    mpi::reduce_Tree_noCoeff(tree, mpi::comm_wrk);
    mpi::broadcast_Tree_noCoeff(tree, mpi::comm_wrk);
}

/** @brief Distribute rank zero function to all ranks */
void mpi::broadcast_function(ComplexFunction &func, MPI_Comm comm) {
/* use same strategy as a reduce, but in reverse order */
#ifdef MRCPP_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) fac *= 2;
    fac /= 2;

    while (fac > 0) {
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 1) {
            // receive
            int src = comm_rank - fac;
            int tag = 4334 + comm_rank;
            mpi::recv_function(func, src, tag, comm);
        }
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 0) {
            // send
            int dst = comm_rank + fac;
            int tag = 4334 + dst;
            if (dst < comm_size) mpi::send_function(func, dst, tag, comm);
        }
        fac /= 2;
    }
    MPI_Barrier(comm);
#endif
}

/** @brief Distribute rank zero function to all ranks */
void mpi::broadcast_Tree_noCoeff(mrcpp::FunctionTree<3> &tree, MPI_Comm comm) {
/* use same strategy as a reduce, but in reverse order */
#ifdef MRCPP_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) fac *= 2;
    fac /= 2;

    while (fac > 0) {
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 1) {
            // receive
            int src = comm_rank - fac;
            int tag = 4334 + comm_rank;
            mrcpp::recv_tree(tree, src, tag, comm, -1, false);
        }
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 0) {
            // send
            int dst = comm_rank + fac;
            int tag = 4334 + dst;
            if (dst < comm_size) mrcpp::send_tree(tree, dst, tag, comm, -1, false);
        }
        fac /= 2;
    }
    MPI_Barrier(comm);
#endif
}

} // namespace mrcpp
