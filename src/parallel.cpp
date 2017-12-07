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
void mpi::communicator::send_tree(FunctionTree<D> &tree, int dst, int tag) {
    SerialFunctionTree<D> &sTree = *tree.getSerialFunctionTree();
    if (sTree.nGenNodes != 0) MSG_FATAL("Sending of GenNodes not implemented");
    
    int n_chunks = sTree.nodeChunks.size();
    MPI_Send(&n_chunks, sizeof(int), MPI_BYTE, dst, tag-1, this->comm);
    println(10, " Sending " << n_chunks << " chunks");

    Timer t1;
    int count = 1;
    for (int iChunk = 0; iChunk < n_chunks; iChunk++) {
        count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
        //set serialIx of the unused nodes to -1
        int iShift = iChunk*sTree.maxNodesPerChunk;
        for (int i = 0; i < sTree.maxNodesPerChunk; i++) {
            if (sTree.nodeStackStatus[iShift+i] !=1)
                sTree.nodeChunks[iChunk][i].setSerialIx(-1);
        }
        MPI_Send(sTree.nodeChunks[iChunk], count, MPI_BYTE, dst, tag+iChunk, this->comm);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        MPI_Send(sTree.nodeCoeffChunks[iChunk], count, MPI_DOUBLE, dst, tag+iChunk+1000, this->comm);
    }
    t1.stop();
    println(10, " Time send                   " << setw(30) << t1);
}
#endif

#ifdef HAVE_MPI
template<int D>
void mpi::communicator::recv_tree(FunctionTree<D> &tree, int src, int tag) {
    MPI_Status status;
    SerialFunctionTree<D> &sTree = *tree.getSerialFunctionTree();

    int n_chunks;
    MPI_Recv(&n_chunks, sizeof(int), MPI_BYTE, src, tag-1, this->comm, &status);
    println(10, " Recieving " << n_chunks << " chunks");

    Timer t1;
    int count = 1;
    for (int iChunk = 0; iChunk < n_chunks; iChunk++) {
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
        MPI_Recv(sTree.nodeChunks[iChunk], count, MPI_BYTE, src, tag+iChunk, this->comm, &status);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        MPI_Recv(sTree.nodeCoeffChunks[iChunk], count, MPI_DOUBLE, src, tag+iChunk+1000, this->comm, &status);
    }
    t1.stop();
    println(10, " Time recieve                " << setw(30) << t1);

    Timer t2;
    sTree.rewritePointers(n_chunks);
    t2.stop();
    println(10, " Time rewrite pointers       " << setw(30) << t2);
}
#endif

#ifdef HAVE_MPI
template<int D>
void mpi::communicator::isend_tree(FunctionTree<D> &tree, int dst, int tag, mpi::requests &reqs) {
    MPI_Request req = MPI_REQUEST_NULL;
    SerialFunctionTree<D> &sTree = *tree.getSerialFunctionTree();
    if (sTree.nGenNodes != 0) MSG_FATAL("Sending of GenNodes not implemented");

    int n_chunks = sTree.nodeChunks.size();
    MPI_Isend(&n_chunks, sizeof(int), MPI_BYTE, dst, tag-1, this->comm, &req);
    println(0, " Sending " << n_chunks << " chunks");

    Timer t1;
    int count = 1;
    for (int iChunk = 0; iChunk < n_chunks; iChunk++) {
        count = sTree.maxNodesPerChunk*sizeof(ProjectedNode<D>);
        //set serialIx of the unused nodes to -1
        int iShift = iChunk*sTree.maxNodesPerChunk;
        for (int i = 0; i < sTree.maxNodesPerChunk; i++) {
            if (sTree.nodeStackStatus[iShift+i] !=1)
                sTree.nodeChunks[iChunk][i].setSerialIx(-1);
        }
        MPI_Isend(sTree.nodeChunks[iChunk], count, MPI_BYTE, dst, tag+iChunk, this->comm, &req);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        MPI_Isend(sTree.nodeCoeffChunks[iChunk], count, MPI_DOUBLE, dst, tag+iChunk+1000, this->comm, &req);
    }
    // Should we push back on individual isends?
    reqs.push_back(req);
    t1.stop();
    println(10, " Time send                   " << setw(30) << t1);
}
#endif

#ifdef HAVE_MPI
void mpi::communicator::wait(mpi::requests &reqs) {
    MPI_Status status;
    Timer t1;
    for (int i = 0; i < reqs.size(); i++) {
        MPI_Wait(&reqs[i], &status);
    }
    t1.stop();
    println(10, " N requests                  " << setw(30) << reqs.size());
    println(10, " Time wait                   " << setw(30) << t1);
    reqs.clear();
}
#endif

#ifdef HAVE_MPI
template void mpi::communicator::send_tree(FunctionTree<1> &tree, int dst, int tag);
template void mpi::communicator::send_tree(FunctionTree<2> &tree, int dst, int tag);
template void mpi::communicator::send_tree(FunctionTree<3> &tree, int dst, int tag);
template void mpi::communicator::recv_tree(FunctionTree<1> &tree, int src, int tag);
template void mpi::communicator::recv_tree(FunctionTree<2> &tree, int src, int tag);
template void mpi::communicator::recv_tree(FunctionTree<3> &tree, int src, int tag);
template void mpi::communicator::isend_tree(FunctionTree<1> &tree, int dst, int tag, mpi::requests &req);
template void mpi::communicator::isend_tree(FunctionTree<2> &tree, int dst, int tag, mpi::requests &req);
template void mpi::communicator::isend_tree(FunctionTree<3> &tree, int dst, int tag, mpi::requests &req);
#endif

