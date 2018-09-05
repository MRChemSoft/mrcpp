#include "mpi_utils.h"
#include "trees/FunctionTree.h"
#include "trees/ProjectedNode.h"
#include "trees/SerialFunctionTree.h"
#include "Printer.h"
#include "Timer.h"

using namespace std;

namespace mrcpp {

/** sh_size in MB
 */
SharedMemory::SharedMemory(MPI_Comm comm, int sh_size) {
#ifdef HAVE_MPI
    int rank;
    MPI_Comm_rank(comm, &rank);
    // MPI_Aint types are used for adresses (can be larger than int)
    MPI_Aint size = (rank == 0) ? 1024*1024*sh_size : 0; //rank 0 defines length of segment
    int disp_unit = 16; //in order for the compiler to keep aligned
    //size is in bytes
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, comm, &this->sh_start_ptr, &this->sh_win);
    MPI_Win_fence(0, this->sh_win); //wait until finished
    MPI_Aint qsize = 0;
    int qdisp = 0;
    MPI_Win_shared_query(this->sh_win, 0, &qsize, &qdisp, &this->sh_start_ptr);
    MPI_Win_fence(0, this->sh_win);
    this->sh_max_ptr = this->sh_start_ptr + qsize/sizeof(double);
    this->sh_end_ptr = this->sh_start_ptr;
#endif
}

SharedMemory::~SharedMemory() {
#ifdef HAVE_MPI
    //deallocates the memory block
    MPI_Win_free(&this->sh_win);
#endif
}

template<int D>
void send_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm comm, int nChunks) {
#ifdef HAVE_MPI
    SerialFunctionTree<D> &sTree = *tree.getSerialFunctionTree();
    if (sTree.nGenNodes != 0) MSG_FATAL("Sending of GenNodes not implemented");

    if (nChunks < 0) {
        nChunks = sTree.getNChunksUsed();
        MPI_Send(&nChunks, sizeof(int), MPI_BYTE, dst, tag, comm);
        println(10, " Sending " << nChunks << " chunks");
    }

    Timer t1;
    int count = 1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
        MPI_Send(sTree.nodeChunks[iChunk], count, MPI_BYTE, dst, tag+iChunk+1, comm);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        MPI_Send(sTree.nodeCoeffChunks[iChunk], count, MPI_DOUBLE, dst, tag+iChunk+1001, comm);
    }
    t1.stop();
    println(10, " Time send                   " << setw(30) << t1);
#endif
}

template<int D>
void recv_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm comm, int nChunks) {
#ifdef HAVE_MPI
    MPI_Status status;
    SerialFunctionTree<D> &sTree = *tree.getSerialFunctionTree();

    if (nChunks < 0) {
        MPI_Recv(&nChunks, sizeof(int), MPI_BYTE, src, tag, comm, &status);
        println(10, " Receiving " << nChunks << " chunks");
    }

    Timer t1;
    int count = 1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        if (iChunk < sTree.nodeChunks.size()) {
            sTree.sNodes = sTree.nodeChunks[iChunk];
        } else {
            double *sNodesCoeff;
            sNodesCoeff = new double[sTree.sizeNodeCoeff*sTree.maxNodesPerChunk];
            sTree.nodeCoeffChunks.push_back(sNodesCoeff);
            sTree.sNodes = (ProjectedNode<D>*) new char[sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>)];
            sTree.nodeChunks.push_back(sTree.sNodes);
        }
        count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
        MPI_Recv(sTree.nodeChunks[iChunk], count, MPI_BYTE, src, tag+iChunk+1, comm, &status);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        MPI_Recv(sTree.nodeCoeffChunks[iChunk], count, MPI_DOUBLE, src, tag+iChunk+1001, comm, &status);
    }
    t1.stop();
    println(10, " Time recieve                " << setw(30) << t1);

    Timer t2;
    sTree.rewritePointers(nChunks);
    t2.stop();
    println(10, " Time rewrite pointers       " << setw(30) << t2);
#endif
}

template<int D>
void isend_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm comm, MPI_Request *req, int nChunks) {
#ifdef HAVE_MPI
    SerialFunctionTree<D> &sTree = *tree.getSerialFunctionTree();
    if (sTree.nGenNodes != 0) MSG_FATAL("Sending of GenNodes not implemented");

    if (nChunks < 0) MSG_FATAL("Number of chunks must be communicated before isend");

    Timer t1;
    int count = 1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
        MPI_Isend(sTree.nodeChunks[iChunk], count, MPI_BYTE, dst, tag+iChunk+1, comm, req);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        MPI_Isend(sTree.nodeCoeffChunks[iChunk], count, MPI_DOUBLE, dst, tag+iChunk+1001, comm, req);
    }
    t1.stop();
    println(10, " Time send                   " << setw(30) << t1);
#endif
}

template<int D>
void share_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm comm) {
#ifdef HAVE_MPI
    Timer t1;
    SerialFunctionTree<D> &sTree = *tree.getSerialFunctionTree();
    if (sTree.nGenNodes != 0) MSG_FATAL("Sending of GenNodes not implemented");

    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    for (int dst = 0; dst < size; dst++) {
        if (dst == src) continue;
        int dst_tag = tag*(dst+1);
        if (rank == src) {
            int nChunks = sTree.nodeChunks.size();
            println(10, " Sending " << nChunks << " chunks");
            MPI_Send(&nChunks, sizeof(int), MPI_BYTE, dst, dst_tag, comm);
            int count = 1;
            for (int iChunk = 0; iChunk < nChunks; iChunk++) {
                count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
                println(10, " Sending chunk " << iChunk);
                MPI_Send(sTree.nodeChunks[iChunk], count, MPI_BYTE, dst, dst_tag+iChunk+1, comm);
                count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
            }
        }
        if (rank == dst) {
            MPI_Status status;

            int nChunks;
            MPI_Recv(&nChunks, sizeof(int), MPI_BYTE, src, dst_tag, comm, &status);
            println(10, " Received " << nChunks << " chunks");

            int count = 1;
            for (int iChunk = 0; iChunk < nChunks; iChunk++) {
                if (iChunk < sTree.nodeChunks.size()) {
                    sTree.sNodes = sTree.nodeChunks[iChunk];
                } else {
                    if (iChunk == 0) sTree.getMemory()->sh_end_ptr = sTree.getMemory()->sh_start_ptr;
                    double *sNodesCoeff = sTree.getMemory()->sh_end_ptr;
                    sTree.getMemory()->sh_end_ptr += (sTree.sizeNodeCoeff*sTree.maxNodesPerChunk);
                    sTree.nodeCoeffChunks.push_back(sNodesCoeff);
                    sTree.sNodes = (ProjectedNode<D>*) new char[sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>)];
                    sTree.nodeChunks.push_back(sTree.sNodes);
                }
                count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
                println(10, " Receiving chunk " << iChunk);
                MPI_Recv(sTree.nodeChunks[iChunk], count, MPI_BYTE, src, dst_tag+iChunk+1, comm, &status);
            }
            sTree.rewritePointers(nChunks);
        }
    }

    t1.stop();
    println(10, " Time share                  " << setw(30) << t1);
#endif
}

template void send_tree(FunctionTree<1> &tree, int dst, int tag, MPI_Comm comm, int nChunks);
template void send_tree(FunctionTree<2> &tree, int dst, int tag, MPI_Comm comm, int nChunks);
template void send_tree(FunctionTree<3> &tree, int dst, int tag, MPI_Comm comm, int nChunks);
template void recv_tree(FunctionTree<1> &tree, int src, int tag, MPI_Comm comm, int nChunks);
template void recv_tree(FunctionTree<2> &tree, int src, int tag, MPI_Comm comm, int nChunks);
template void recv_tree(FunctionTree<3> &tree, int src, int tag, MPI_Comm comm, int nChunks);
template void isend_tree(FunctionTree<1> &tree, int dst, int tag, MPI_Comm comm, MPI_Request *req, int nChunks);
template void isend_tree(FunctionTree<2> &tree, int dst, int tag, MPI_Comm comm, MPI_Request *req, int nChunks);
template void isend_tree(FunctionTree<3> &tree, int dst, int tag, MPI_Comm comm, MPI_Request *req, int nChunks);
template void share_tree(FunctionTree<1> &tree, int src, int tag, MPI_Comm comm);
template void share_tree(FunctionTree<2> &tree, int src, int tag, MPI_Comm comm);
template void share_tree(FunctionTree<3> &tree, int src, int tag, MPI_Comm comm);

} // namespace mrcpp
