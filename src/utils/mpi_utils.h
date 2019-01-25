#pragma once

#ifdef HAVE_MPI
#include <mpi.h>
#else
using MPI_Comm = int;
using MPI_Win = int;
using MPI_Request = int;
#endif

namespace mrcpp {

/** Share memory within a compute node
 */
class SharedMemory {
public:
    SharedMemory(MPI_Comm comm, int sh_size);
    SharedMemory(const SharedMemory &mem) = delete;
    SharedMemory &operator=(const SharedMemory &mem) = delete;
    ~SharedMemory();

    double *sh_start_ptr;  //start of shared block
    double *sh_end_ptr;    //end of used part
    double *sh_max_ptr;    //end of shared block
    MPI_Win sh_win;        //MPI window object
    int rank;              //rank among shared group
};

template<int D> class FunctionTree;

template<int D> void isend_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm comm, MPI_Request *req, int nChunks = -1);
template<int D> void send_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm comm, int nChunks = -1);
template<int D> void recv_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm comm, int nChunks = -1);
template<int D> void share_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm comm);

}

