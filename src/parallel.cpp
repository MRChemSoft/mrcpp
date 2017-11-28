#include "parallel.h"
#include "FunctionTree.h"
#include "ProjectedNode.h"
#include "SerialFunctionTree.h"
#include "Printer.h"
#include "Timer.h"

mpi::communicator mpi::world;

using namespace std;

#ifdef HAVE_MPI
template<int D>
void mpi::communicator::send_tree(mpi::tree_request<D> &request, int dst) {
    if (request.s_tree == 0) MSG_FATAL("Invalid request");

    SerialFunctionTree<D> &sTree = *request.s_tree;
    if (sTree.nGenNodes != 0) MSG_FATAL("Sending of GenNodes not implemented");
    
    MPI_Send(&request.n_chunks, sizeof(int), MPI_BYTE, dst, request.tag, this->comm);
    println(10, " Sending " << request.n_chunks << " chunks");

    Timer t1;
    int count = 1;
    for (int iChunk = 0; iChunk < request.n_chunks; iChunk++) {
        count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
        //set serialIx of the unused nodes to -1
        int iShift = iChunk*sTree.maxNodesPerChunk;
        for (int i = 0; i < sTree.maxNodesPerChunk; i++) {
            if (sTree.nodeStackStatus[iShift+i] !=1)
                sTree.nodeChunks[iChunk][i].setSerialIx(-1);
        }
        MPI_Send(sTree.nodeChunks[iChunk], count, MPI_BYTE, dst, request.tag+1+iChunk, this->comm);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        MPI_Send(sTree.nodeCoeffChunks[iChunk], count, MPI_DOUBLE, dst, request.tag+iChunk*1000, this->comm);
    }
    t1.stop();
    println(10, " Time send                   " << setw(30) << t1);
}
#endif

#ifdef HAVE_MPI
template<int D>
void mpi::communicator::recv_tree(mpi::tree_request<D> &request, int src) {
    if (request.s_tree == 0) MSG_FATAL("Invalid request");

    SerialFunctionTree<D> &sTree = *request.s_tree;

    MPI_Recv(&request.n_chunks, sizeof(int), MPI_BYTE, src, request.tag, this->comm, &request.mpi_stat);
    println(10, " Recieving " << request.n_chunks << " chunks");

    Timer t1;
    int count = 1;
    for (int iChunk = 0; iChunk < request.n_chunks; iChunk++) {
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
        MPI_Recv(sTree.nodeChunks[iChunk], count, MPI_BYTE, src, request.tag+1+iChunk, this->comm, &request.mpi_stat);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        MPI_Recv(sTree.nodeCoeffChunks[iChunk], count, MPI_DOUBLE, src, request.tag+iChunk*1000, this->comm, &request.mpi_stat);
    }
    t1.stop();
    println(10, " Time recieve                " << setw(30) << t1);

    Timer t2;
    sTree.rewritePointers(request.n_chunks);
    t2.stop();
    println(10, " Time rewrite pointers       " << setw(30) << t2);
}
#endif

#ifdef HAVE_MPI
template void mpi::communicator::send_tree(mpi::tree_request<1> &request, int dst);
template void mpi::communicator::send_tree(mpi::tree_request<2> &request, int dst);
template void mpi::communicator::send_tree(mpi::tree_request<3> &request, int dst);
template void mpi::communicator::recv_tree(mpi::tree_request<1> &request, int src);
template void mpi::communicator::recv_tree(mpi::tree_request<2> &request, int src);
template void mpi::communicator::recv_tree(mpi::tree_request<3> &request, int src);
#endif

