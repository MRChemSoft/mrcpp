#pragma once

#ifdef HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
typedef int MPI_Win;
typedef int MPI_Request;
#endif

/** Share memory within a compute node
 */
class SharedMemory {
public:
    SharedMemory(MPI_Comm comm, int sh_size);
    ~SharedMemory();

    double *sh_start_ptr;  //start of shared block
    double *sh_end_ptr;    //end of used part
    double *sh_max_ptr;    //end of shared block
    MPI_Win sh_win;        //MPI window object
};

template<int D> class FunctionTree;

namespace mrcpp {
template<int D> void isend_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm comm, MPI_Request *req);
template<int D> void send_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm comm);
template<int D> void recv_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm comm);
template<int D> void share_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm comm);
}

