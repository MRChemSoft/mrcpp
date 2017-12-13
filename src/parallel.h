#pragma once

#include <vector>
#include "config.h"
#include "constants.h"

#define EIGEN_DONT_PARALLELIZE

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_set_dynamic(n)
#define omp_set_lock(x)
#define omp_unset_lock(x)
#define omp_test_lock(x)
#endif

template<int D> class FunctionTree;

#ifdef HAVE_MPI

/** Share memory within a compute node
 */
class SharedMemory {
public:
    SharedMemory(MPI_Comm &comm, int sh_size = 500);
    ~SharedMemory();

    double *sh_start_ptr;  //start of shared block
    double *sh_end_ptr;    //end of used part
    double *sh_max_ptr;    //end of shared block
    MPI_Win sh_win;        //MPI window object
};

namespace mpi {

template<int D> void isend_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm &comm, MPI_Request &req);
template<int D> void send_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm &comm);
template<int D> void recv_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm &comm);
template<int D> void share_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm &comm);

};

#else

typedef int MPI_Comm;

class SharedMemory {
public:
    SharedMemory(MPI_Comm &comm, int sh_size = 0) { }
    ~SharedMemory() { }
};

namespace mpi {

template<int D> void isend_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm &comm, MPI_Request &req) { }
template<int D> void send_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm &comm) { }
template<int D> void recv_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm &comm) { }
template<int D> void share_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm &comm) { }

};

#endif

