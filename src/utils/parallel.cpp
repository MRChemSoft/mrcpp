/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include <Eigen/Core>
#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <thread>
#include <vector>

#include "Bank.h"
#include "omp_utils.h"
#include "parallel.h"
#include "trees/FunctionTree.h"

#ifdef MRCPP_HAS_OMP
#define mrcpp_get_max_threads() omp_get_max_threads()
#define mrcpp_get_num_procs() omp_get_num_procs()
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
int is_bankmaster = 0;
int bank_size = 0;
int bank_per_node = 0;
int omp_threads = -1;
int use_omp_num_threads = -1;
int tot_bank_size = 0;
int max_tag = 0;
vector<int> bankmaster;
int task_bank = -1;

MPI_Comm comm_wrk;
MPI_Comm comm_share;
MPI_Comm comm_sh_group;
MPI_Comm comm_bank;

int id_shift;

extern int metadata_block[3];
extern int const size_metadata = 3;

void initialize() {
    Eigen::setNbThreads(1);
    mrcpp_set_dynamic(0);

#ifdef MRCPP_HAS_MPI
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Comm node_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
    int node_rank, node_size;
    MPI_Comm_rank(node_comm, &node_rank);
    MPI_Comm_size(node_comm, &node_size);

    comm_bank = MPI_COMM_WORLD;
    MPI_Comm comm_remainder;

    if (world_size < 2) {
        bank_size = 0;
    } else if (bank_size < 0) {
        if (bank_per_node >= 0) {
            bank_size = node_size * bank_per_node;
        } else {
            bank_size = max(world_size / 3, 1);
        }
    } else if (bank_size >= 0 and bank_per_node >= 0) {
        if (bank_size != node_size * bank_per_node and world_rank == 0) std::cout << "WARNING: bank_size and bank_per_node are incompatible " << bank_size << " " << bank_per_node << std::endl;
    }
    if (world_size - bank_size < 1) MSG_ABORT("No MPI ranks left for working!");
    if (bank_size < 1 and world_size > 1) MSG_ABORT("Bank size must be at least one when using MPI!");

    bankmaster.resize(bank_size);
    for (int i = 0; i < bank_size; i++) {
        bankmaster[i] = world_size - i - 1;
    }
    if (world_rank < world_size - bank_size) {
        is_bank = 0;
        is_centralbank = 0;
        is_bankclient = 1;
    } else {
        is_bank = 1;
        is_centralbank = 1;
        is_bankclient = 0;
        if (world_rank == world_size - bank_size) is_bankmaster = 1;
    }
    MPI_Comm_split(MPI_COMM_WORLD, is_bankclient, world_rank, &comm_remainder);

    MPI_Comm_split_type(comm_remainder, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_share);

    MPI_Comm_rank(comm_share, &share_rank);
    MPI_Comm_size(comm_share, &share_size);

    MPI_Comm_split(comm_remainder, share_rank, world_rank, &comm_sh_group);

    MPI_Comm_rank(comm_sh_group, &sh_group_rank);

    wrk_rank = share_rank + sh_group_rank * world_size;
    MPI_Comm_split(comm_remainder, 0, wrk_rank, &comm_wrk);

    MPI_Comm_rank(comm_wrk, &wrk_rank);
    MPI_Comm_size(comm_wrk, &wrk_size);

    tot_bank_size = bank_size;
    if (bank_size <= 2 and bank_size > 0) {
        task_bank = bankmaster[0];
    } else if (bank_size > 1) {
        bank_size--;
        task_bank = bankmaster[bank_size];
    }

    void *val;
    int flag;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &val, &flag);
    max_tag = *(int *)val / 2;
    id_shift = max_tag / 2;

    MPI_Comm comm_share_world;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_share_world);

    int n_bank_thisnode;
    MPI_Allreduce(&is_bank, &n_bank_thisnode, 1, MPI_INT, MPI_SUM, comm_share_world);
    int n_wrk_thisnode;
    MPI_Allreduce(&is_bankclient, &n_wrk_thisnode, 1, MPI_INT, MPI_SUM, comm_share_world);

    int omp_threads_available = thread::hardware_concurrency();

    int nthreads = 1;
    int my_OMP_NUM_THREADS = mrcpp_get_max_threads();
    MPI_Bcast(&my_OMP_NUM_THREADS, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (use_omp_num_threads) {
        int total_omp_threads_per_node = my_OMP_NUM_THREADS * (n_bank_thisnode + n_wrk_thisnode);
        nthreads = (total_omp_threads_per_node - n_bank_thisnode) / n_wrk_thisnode;
    } else {
        if (is_bankclient) nthreads = (omp_threads_available / 2 - n_bank_thisnode) / n_wrk_thisnode;
        nthreads = min(nthreads, mrcpp_get_num_procs() / 2);
        if (is_bank) nthreads = 1;
        if (omp_threads > 0) {
            if (omp_threads != nthreads and world_rank == 0) {
                cout << "Warning: recommended number of threads is " << nthreads << endl;
                cout << "setting number of threads to omp_threads, " << max(1, omp_threads) << endl;
            }
            nthreads = omp_threads;
        }
    }
    nthreads = max(1, nthreads);

    if (nthreads * n_wrk_thisnode + n_bank_thisnode < omp_threads_available / 3 and world_rank == 0) {
        std::cout << "WARNING: only " << nthreads * n_wrk_thisnode + n_bank_thisnode << " threads used per node while " << omp_threads_available << " logical cpus are accessible " << std::endl;
    }

    if (nthreads > mrcpp_get_num_procs() / 2) { std::cout << "WARNING: MPI rank " << world_rank << " will use " << nthreads << " but only " << mrcpp_get_num_procs() / 2 << " procs are accessible" << std::endl; }

    omp::n_threads = nthreads;
    mrcpp::set_max_threads(nthreads);

    if (is_bank) {
        if (is_centralbank) { dataBank.open(); }
        finalize();
        exit(EXIT_SUCCESS);
    }
#else
    bank_size = 0;
    mrcpp::set_max_threads(omp::n_threads);
#endif
}

void finalize() {
#ifdef MRCPP_HAS_MPI
    if (bank_size > 0 and grand_master()) {
        println(4, " max data in bank " << dataBank.get_maxtotalsize() << " MB ");
        dataBank.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
}

void barrier(MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    MPI_Barrier(comm);
#endif
}

bool grand_master() {
    return (world_rank == 0 and is_bankclient) ? true : false;
}

bool share_master() {
    return (share_rank == 0) ? true : false;
}

bool my_func(int j) {
    return ((j) % wrk_size == wrk_rank) ? true : false;
}

bool my_func(const CompFunction<3> &func) {
    return my_func(func.rank());
}

bool my_func(CompFunction<3> *func) {
    return my_func(func->rank());
}

void free_foreign(CompFunctionVector &Phi) {
    for (CompFunction<3> &i : Phi) {
        if (not my_func(i)) i.free();
    }
}

void allreduce_vector(IntVector &vec, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

void allreduce_vector(DoubleVector &vec, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

void allreduce_vector(ComplexVector &vec, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

void allreduce_matrix(IntMatrix &mat, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

void allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

void allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

void send_function(const CompFunction<3> &func, int dst, int tag, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    for (int i = 0; i < func.Ncomp(); i++) {
        if (func.isreal())
            func.Nchunks()[i] = func.CompD[i]->getNChunks();
        else
            func.Nchunks()[i] = func.CompC[i]->getNChunks();
    }
    MPI_Send(&func.func_ptr->data, sizeof(CompFunctionData<3>), MPI_BYTE, dst, 0, comm);
    for (int i = 0; i < func.Ncomp(); i++) {
        if (func.isreal())
            mrcpp::send_tree(*func.CompD[i], dst, tag, comm, func.Nchunks()[i]);
        else
            mrcpp::send_tree(*func.CompC[i], dst, tag, comm, func.Nchunks()[i]);
    }
#endif
}

void recv_function(CompFunction<3> &func, int src, int tag, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int func_ncomp_in = func.Ncomp();
    MPI_Recv(&func.func_ptr->data, sizeof(CompFunctionData<3>), MPI_BYTE, src, 0, comm, &status);
    for (int i = 0; i < func.Ncomp(); i++) {
        if (func_ncomp_in <= i) func.alloc(i + 1);
        if (func.isreal())
            mrcpp::recv_tree(*func.CompD[i], src, tag, comm, func.Nchunks()[i]);
        else
            mrcpp::recv_tree(*func.CompC[i], src, tag, comm, func.Nchunks()[i]);
    }
#endif
}

void share_function(CompFunction<3> &func, int src, int tag, MPI_Comm comm) {
    if (func.isShared()) {
#ifdef MRCPP_HAS_MPI
        for (int comp = 0; comp < func.Ncomp(); comp++) {
            if (func.isreal())
                mrcpp::share_tree(*func.CompD[comp], src, tag, comm);
            else
                mrcpp::share_tree(*func.CompC[comp], src, tag, comm);
        }
#endif
    }
}

void reduce_function(double prec, CompFunction<3> &func, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1;
    while (fac < comm_size) {
        if ((comm_rank / fac) % 2 == 0) {
            int src = comm_rank + fac;
            if (src < comm_size) {
                CompFunction<3> func_i;
                int tag = 3333 + src;
                recv_function(func_i, src, tag, comm);
                func.add(1.0, func_i);
                func.crop(prec);
            }
        }
        if ((comm_rank / fac) % 2 == 1) {
            int dest = comm_rank - fac;
            if (dest >= 0) {
                int tag = 3333 + comm_rank;
                send_function(func, dest, tag, comm);
                break;
            }
        }
        fac *= 2;
    }
    MPI_Barrier(comm);
#endif
}

template <typename T> void reduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1;
    while (fac < comm_size) {
        if ((comm_rank / fac) % 2 == 0) {
            int src = comm_rank + fac;
            if (src < comm_size) {
                int tag = 3333 + src;
                mrcpp::FunctionTree<3, T> tree_i(tree.getMRA());
                mrcpp::recv_tree(tree_i, src, tag, comm, -1, false);
                tree.appendTreeNoCoeff(tree_i);
            }
        }
        if ((comm_rank / fac) % 2 == 1) {
            int dest = comm_rank - fac;
            if (dest >= 0) {
                int tag = 3333 + comm_rank;
                mrcpp::send_tree(tree, dest, tag, comm, -1, false);
                break;
            }
        }
        fac *= 2;
    }
    MPI_Barrier(comm);
#endif
}

template <typename T> void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, vector<FunctionTree<3, T>> &Phi, MPI_Comm comm) {
    int N = Phi.size();
    for (int j = 0; j < N; j++) {
        if (not my_func(j)) continue;
        tree.appendTreeNoCoeff(Phi[j]);
    }
#ifdef MRCPP_HAS_MPI
    mrcpp::mpi::reduce_Tree_noCoeff(tree, comm_wrk);
    mrcpp::mpi::broadcast_Tree_noCoeff(tree, comm_wrk);
#endif
}

template <typename T> void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, vector<CompFunction<3>> &Phi, MPI_Comm comm) {
    int N = Phi.size();
    for (int j = 0; j < N; j++) {
        if (not my_func(j)) continue;
        if (Phi[j].isreal()) tree.appendTreeNoCoeff(*Phi[j].CompD[0]);
        if (Phi[j].iscomplex()) tree.appendTreeNoCoeff(*Phi[j].CompC[0]);
    }
#ifdef MRCPP_HAS_MPI
    mrcpp::mpi::reduce_Tree_noCoeff(tree, comm_wrk);
    mrcpp::mpi::broadcast_Tree_noCoeff(tree, comm_wrk);
#endif
}

void broadcast_function(CompFunction<3> &func, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1;
    while (fac < comm_size) fac *= 2;
    fac /= 2;

    while (fac > 0) {
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 1) {
            int src = comm_rank - fac;
            int tag = 4334 + comm_rank;
            recv_function(func, src, tag, comm);
        }
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 0) {
            int dst = comm_rank + fac;
            int tag = 4334 + dst;
            if (dst < comm_size) send_function(func, dst, tag, comm);
        }
        fac /= 2;
    }
    MPI_Barrier(comm);
#endif
}

template <typename T> void broadcast_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1;
    while (fac < comm_size) fac *= 2;
    fac /= 2;

    while (fac > 0) {
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 1) {
            int src = comm_rank - fac;
            int tag = 4334 + comm_rank;
            mrcpp::recv_tree(tree, src, tag, comm, -1, false);
        }
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 0) {
            int dst = comm_rank + fac;
            int tag = 4334 + dst;
            if (dst < comm_size) mrcpp::send_tree(tree, dst, tag, comm, -1, false);
        }
        fac /= 2;
    }
    MPI_Barrier(comm);
#endif
}

template void reduce_Tree_noCoeff(mrcpp::FunctionTree<3, double> &tree, MPI_Comm comm);
template void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, double> &tree, std::vector<FunctionTree<3, double>> &Phi, MPI_Comm comm);
template void broadcast_Tree_noCoeff(mrcpp::FunctionTree<3, double> &tree, MPI_Comm comm);
template void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, double> &tree, std::vector<CompFunction<3>> &Phi, MPI_Comm comm);

template void reduce_Tree_noCoeff(mrcpp::FunctionTree<3, ComplexDouble> &tree, MPI_Comm comm);
template void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, ComplexDouble> &tree, std::vector<FunctionTree<3, ComplexDouble>> &Phi, MPI_Comm comm);
template void broadcast_Tree_noCoeff(mrcpp::FunctionTree<3, ComplexDouble> &tree, MPI_Comm comm);
template void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, ComplexDouble> &tree, std::vector<CompFunction<3>> &Phi, MPI_Comm comm);

} // namespace mpi
} // namespace mrcpp
