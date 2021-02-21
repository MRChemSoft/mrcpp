/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "FunctionNodeAllocator.h"

#include <stack>

#include "FunctionTree.h"
#include "FunctionNode.h"
#include "utils/Printer.h"
#include "utils/mpi_utils.h"

namespace mrcpp {

/** SerialTree class constructor.
 * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
 * Initializes rootNodes to represent the zero function and allocate their nodes.
 * NOTES:
 * Serial trees are made of projected nodes, and include gennodes and loose nodes separately.
 * All created (using class creator) Projected nodes or GenNodes are loose nodes.
 * Loose nodes have their coeff in serial Tree, but not the node part.
 * Projected nodes and GenNodes that are created by their creator, are detroyed by destructor ~FunctionNode and
 * ~GenNode. Serial tree nodes are not using the destructors, but explicitely call to deallocNodes or deallocGenNodes
 * Gen nodes and loose nodes are not counted with MWTree->[in/de]crementNodeCount()
 */
template <int D>
FunctionNodeAllocator<D>::FunctionNodeAllocator(FunctionTree<D> *tree, SharedMemory *mem, bool gen)
        : NodeAllocator<D>(tree, mem)
        , genNode(gen)
        , lastNode(nullptr) {
    // reserve space for chunk pointers to avoid excessive reallocation
    this->nodeChunks.reserve(100);
    this->nodeCoeffChunks.reserve(100);

    int tDim = (this->genNode) ? 1 : (1 << D); // genNodes have only one block
    this->coeffsPerNode = tDim * this->tree_p->getKp1_d();
    println(10, "SizeNode Coeff (kB) " << this->coeffsPerNode * sizeof(double) / 1024);

    if (D < 3) {
        // define from number of nodes per chunk
        this->maxNodesPerChunk = 64;
    } else {
        // 2 MB small for no waisting place, but large enough so that latency and overhead work is negligible
        int sizePerChunk = 2 * 1024 * 1024;
        this->maxNodesPerChunk = (sizePerChunk / this->coeffsPerNode / sizeof(double) / 8) * 8;
    }

    MRCPP_INIT_OMP_LOCK();
}

template <int D> FunctionNodeAllocator<D>::~FunctionNodeAllocator() {
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[](char *)(this->nodeChunks[i]);
    if (not this->isShared()) // if the data is shared, it must be freed by MPI_Win_free
        for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];

    this->nodeStackStatus.clear();
    MRCPP_DESTROY_OMP_LOCK();
}

/** reset the start node counter */
template <int D> void FunctionNodeAllocator<D>::clear(int n) {
    MRCPP_SET_OMP_LOCK();
    std::fill(this->nodeStackStatus.begin() + n, this->nodeStackStatus.end(), 0);
    this->topStack = n;
    this->lastNode = getNodeNoLock(n);

    if (this->isShared()) {
        this->shmem_p->sh_end_ptr = this->shmem_p->sh_start_ptr + (n / this->maxNodesPerChunk + 1) * this->coeffsPerNode * this->maxNodesPerChunk;
    }
    MRCPP_UNSET_OMP_LOCK();
}

template <int D> FunctionNode<D> * FunctionNodeAllocator<D>::getNode_p(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int cIdx = sIdx % this->maxNodesPerChunk;   // position in chunk
    FunctionNode<D> *node = this->nodeChunks[chunk] + cIdx;
    MRCPP_UNSET_OMP_LOCK();
    return node;
}

template <int D> double * FunctionNodeAllocator<D>::getCoef_p(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int idx = sIdx % this->maxNodesPerChunk;   // position in chunk
    double *coefs = this->nodeCoeffChunks[chunk] + idx * this->coeffsPerNode;
    MRCPP_UNSET_OMP_LOCK();
    return coefs;
}

template <int D> FunctionNode<D> * FunctionNodeAllocator<D>::getNodeNoLock(int sIdx) {
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int cIdx = sIdx % this->maxNodesPerChunk;   // position in chunk
    FunctionNode<D> *node = this->nodeChunks[chunk] + cIdx;
    return node;
}

template <int D> double * FunctionNodeAllocator<D>::getCoefNoLock(int sIdx) {
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int idx = sIdx % this->maxNodesPerChunk;   // position in chunk
    double *coefs = this->nodeCoeffChunks[chunk] + idx * this->coeffsPerNode;
    return coefs;
}

template <int D> int FunctionNodeAllocator<D>::alloc(int nAlloc, bool coeff) {
    MRCPP_SET_OMP_LOCK();

    // move topstack to start of next chunk if current chunk is too small
    int cIdx = this->topStack % (this->maxNodesPerChunk);
    bool chunkOverflow = ((cIdx + nAlloc) > this->maxNodesPerChunk);
    if (chunkOverflow) this->topStack = this->maxNodesPerChunk * ((this->topStack + nAlloc - 1) / this->maxNodesPerChunk);

    // append chunk if necessary
    int chunk = this->topStack / this->maxNodesPerChunk;
    bool needNewChunk = (chunk >= this->nodeChunks.size());
    if (needNewChunk) appendChunk(coeff);

    // return value is index of first new node
    auto sIdx = this->topStack;

    // fill stack status
    auto &status = this->nodeStackStatus;
    for (int i = sIdx; i < sIdx + nAlloc; i++) {
        if (status[i] != 0) MSG_ERROR(" NodeStackStatus: not available [" << i << "] : " << status[i]);
        status[i] = 1;
    }

    // advance stack pointers
    this->nNodes += nAlloc;
    this->topStack += nAlloc;
    this->lastNode = getNodeNoLock(sIdx) + nAlloc;
    MRCPP_UNSET_OMP_LOCK();

    return sIdx;
}

template <int D> void FunctionNodeAllocator<D>::dealloc(int serialIx) {
    MRCPP_SET_OMP_LOCK();
    this->nodeStackStatus[serialIx] = 0; // mark as available
    if (serialIx == this->topStack - 1) {  // top of stack
        while (this->nodeStackStatus[this->topStack - 1] == 0) {
            this->topStack--;
            if (this->topStack < 1) break;
        }
        // has to redefine lastNode
        this->lastNode = getNodeNoLock(this->topStack);
    }
    this->nNodes--;
    MRCPP_UNSET_OMP_LOCK();
}

template <int D> void FunctionNodeAllocator<D>::initChunk(int iChunk, bool coeff) {
    MRCPP_SET_OMP_LOCK();
    if (iChunk >= getNChunks()) appendChunk(coeff);
    MRCPP_UNSET_OMP_LOCK();
}

template <int D> void FunctionNodeAllocator<D>::appendChunk(bool coeff) {
    // make coeff chunk
    if (coeff) {
        double *c_chunk = nullptr;
        if (this->isShared()) {
            // for coefficients, take from the shared memory block
            c_chunk = this->shmem_p->sh_end_ptr;
            this->shmem_p->sh_end_ptr += (this->coeffsPerNode * this->maxNodesPerChunk);
            // may increase size dynamically in the future
            if (this->shmem_p->sh_max_ptr < this->shmem_p->sh_end_ptr) MSG_ABORT("Shared block too small");
        } else {
            c_chunk = new double[getCoeffChunkSize()];
        }
        this->nodeCoeffChunks.push_back(c_chunk);
    }

    // make node chunk
    auto n_chunk = (FunctionNode<D> *)new char[getNodeChunkSize()];
    for (int i = 0; i < this->maxNodesPerChunk; i++) {
        n_chunk[i].serialIx = -1;
        n_chunk[i].parentSerialIx = -1;
        n_chunk[i].childSerialIx = -1;
    }
    this->nodeChunks.push_back(n_chunk);

    // append to nodeStackStatus
    int oldsize = this->nodeStackStatus.size();
    int newsize = oldsize + this->maxNodesPerChunk;
    this->nodeStackStatus.resize(newsize);
    std::fill(this->nodeStackStatus.begin() + oldsize, this->nodeStackStatus.end(), 0);
}

/** Fill all holes in the chunks with occupied nodes, then remove all empty chunks */
template <int D> int FunctionNodeAllocator<D>::shrinkChunks() {
    MRCPP_SET_OMP_LOCK();
    int nAlloc = (1 << D);
    if (this->maxNodesPerChunk * this->nodeChunks.size() <=
        this->getTree()->getNNodes() + this->maxNodesPerChunk + nAlloc - 1) {
        MRCPP_UNSET_OMP_LOCK();
        return 0; // no chunks to remove
    }

    int posocc = 0;
    int posavail = this->getTree()->getRootBox().size(); // start after root nodes
    while (true) {
        posavail = findNextAvailable(posavail, nAlloc);
        if (posavail >= this->topStack) break; // treated all nodes

        posocc = findNextOccupied(posavail);
        if (posocc >= this->topStack) break; // treated all nodes

        moveNodes(nAlloc, posocc, posavail);
    }

    // find the last used node
    posocc = this->topStack - 1;
    while (this->nodeStackStatus[posocc] == 0 and posocc > 0) posocc--;
    this->topStack = posocc + 1;
    this->lastNode = getNodeNoLock(this->topStack);

    int nChunksDeleted = deleteUnusedChunks();
    this->getTree()->resetEndNodeTable();

    MRCPP_UNSET_OMP_LOCK();
    return nChunksDeleted;
}

template <int D> int FunctionNodeAllocator<D>::findNextAvailable(int pos, int nAlloc) const {
    // Last positions on a chunk cannot be used if there is no place for nAlloc siblings on the same chunk
    bool chunkTooSmall = (pos + nAlloc - 1) / this->maxNodesPerChunk != pos / this->maxNodesPerChunk;
    bool posIsOccupied = (this->nodeStackStatus[pos] != 0);
    bool endOfStack = (pos >= this->topStack);
    while ((posIsOccupied or chunkTooSmall) and not(endOfStack)) {
        pos++;
        chunkTooSmall = (pos + nAlloc - 1) / this->maxNodesPerChunk != pos / this->maxNodesPerChunk;
        posIsOccupied = (this->nodeStackStatus[pos] != 0);
        endOfStack = (pos >= this->topStack);
    }
    return pos;
}

template <int D> int FunctionNodeAllocator<D>::findNextOccupied(int pos) const {
    bool posIsAvailable = (this->nodeStackStatus[pos] == 0);
    bool endOfStack = (pos >= this->topStack);
    while (posIsAvailable and not(endOfStack)) {
        pos++;
        posIsAvailable = (this->nodeStackStatus[pos] == 0);
        endOfStack = (pos >= this->topStack);
    }
    return pos;
}

template <int D> void FunctionNodeAllocator<D>::moveNodes(int nNodes, int srcIdx, int dstIdx) {
    FunctionNode<D> *srcNode = getNodeNoLock(srcIdx);
    FunctionNode<D> *dstNode = getNodeNoLock(dstIdx);

    // check that all siblings are consecutive. Should never be root node.
    for (int i = 0; i < nNodes; i++) assert(this->nodeStackStatus[dstIdx + i] == 0);
    for (int i = 1; i < nNodes; i++) assert((srcNode + i)->parent->serialIx == srcNode->parent->serialIx); // siblings

    // just copy everything "as is"
    for (int i = 0; i < nNodes * sizeof(FunctionNode<D>); i++) ((char *)dstNode)[i] = ((char *)srcNode)[i];

    // coefs have new adresses
    double *coefs_p = getCoefNoLock(dstIdx);
    for (int i = 0; i < nNodes; i++) (dstNode + i)->coefs = coefs_p + i * this->getNCoefs();

    // copy coefs to new adress
    if (not this->isShared()) {
        for (int i = 0; i < nNodes * this->coeffsPerNode; i++) dstNode->coefs[i] = srcNode->coefs[i];
    } else {
        if (this->shmem_p->rank == 0) // only master copy the data. careful with sync
            for (int i = 0; i < nNodes * this->coeffsPerNode; i++) dstNode->coefs[i] = srcNode->coefs[i];
    }

    // update node
    for (int i = 0; i < nNodes; i++) (dstNode + i)->serialIx = dstIdx + i;

    // update parent
    dstNode->parent->childSerialIx = dstIdx;
    for (int i = 0; i < nNodes; i++) dstNode->parent->children[i] = dstNode + i;

    // update children
    for (int i = 0; i < nNodes; i++) {
        for (int j = 0; j < (dstNode + i)->getNChildren(); j++) {
            (dstNode + i)->children[j]->parentSerialIx = dstIdx + i;
            (dstNode + i)->children[j]->parent = dstNode + i;
        }
    }

    // mark moved nodes as occupied
    for (int i = 0; i < nNodes; i++) this->nodeStackStatus[dstIdx + i] = 1;
    dstIdx += nNodes;

    // delete "old" nodes
    for (int i = 0; i < nNodes; i++) this->nodeStackStatus[srcIdx + i] = 0;
    for (int i = 0; i < nNodes; i++) (srcNode + i)->serialIx = -1;
}

template <int D> int FunctionNodeAllocator<D>::deleteUnusedChunks() {
    // number of occupied chunks
    int nChunksTotal = getNChunks();
    int nChunksUsed = getNChunksUsed();
    for (int i = nChunksUsed; i < this->nodeChunks.size(); i++) delete[](char *)(this->nodeChunks[i]);

    if (this->isShared()) {
        // shared coefficients cannot be fully deallocated, only pointer is moved.
        this->shmem_p->sh_end_ptr -= (nChunksTotal - nChunksUsed) * this->coeffsPerNode * this->maxNodesPerChunk;
    } else {
        for (int i = nChunksUsed; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];
    }

    // shrink the stacks
    this->nodeChunks.resize(nChunksUsed);
    this->nodeCoeffChunks.resize(nChunksUsed);
    this->nodeStackStatus.resize(nChunksUsed * this->maxNodesPerChunk);
    return nChunksTotal - nChunksUsed;
}

/** Traverse tree and redefine pointer, counter and tables. */
template <int D> void FunctionNodeAllocator<D>::rewritePointers(bool coeff) {
    MRCPP_SET_OMP_LOCK();
    this->nNodes = 0;
    this->getTree()->nodesAtDepth.clear();
    this->getTree()->squareNorm = 0.0;
    this->getTree()->clearEndNodeTable();

    // make virtual table pointers
    char *cvptr_FunctionNode = nullptr;
    {
        FunctionNode<D> tmp;
        cvptr_FunctionNode = *(char **)(&tmp);
    }

    // reinitialize stacks
    // nodeChunks have been adapted to receiving tree and maybe larger than nodeStackStatus
    int nodeCount = this->nodeChunks.size() * this->maxNodesPerChunk;
    this->nodeStackStatus.resize(nodeCount);
    std::fill(this->nodeStackStatus.begin(), this->nodeStackStatus.end(), 0);

    NodeBox<D> &rootbox = this->getTree()->getRootBox();
    MWNode<D> **roots = rootbox.getNodes();

    std::stack<FunctionNode<D> *> stack;
    for (int rIdx = 0; rIdx < rootbox.size(); rIdx++) {
        auto *root_p = getNodeNoLock(rIdx);
        stack.push(root_p);
        roots[rIdx] = root_p;
    }
    this->topStack = 0;
    while (not stack.empty()) {
        auto *node = stack.top();
        auto sIdx = node->serialIx;
        auto pIdx = node->parentSerialIx;
        auto cIdx = node->childSerialIx;

        this->nNodes++;
        this->topStack = std::max(this->topStack, sIdx + 1);
        this->getTree()->incrementNodeCount(node->getScale());
        if (node->isEndNode()) this->getTree()->squareNorm += node->getSquareNorm();
        if (node->isEndNode()) this->getTree()->endNodeTable.push_back(node);

        // normally (intel) the virtual table does not change, but we overwrite anyway
        *(char **)(node) = cvptr_FunctionNode;

        node->tree = this->getTree();
        node->coefs = (coeff) ? this->getCoefNoLock(sIdx) : nullptr;
        node->parent = (pIdx >= 0) ? this->getNodeNoLock(pIdx) : nullptr;

        stack.pop();
        auto *child_p = this->getNodeNoLock(cIdx);
        for (int i = 0; i < node->getNChildren(); i++) {
            node->children[i] = child_p;
            stack.push(child_p);
            child_p++;
        }
        this->nodeStackStatus[sIdx] = 1; // occupied
    }
    this->lastNode = this->getNodeNoLock(this->topStack);
    MRCPP_UNSET_OMP_LOCK();
}

template <int D> void FunctionNodeAllocator<D>::print() const {
    int n = 0;
    for (int iChunk = 0; iChunk < getNChunks(); iChunk++) {
        int iShift = iChunk * this->maxNodesPerChunk;
        printout(0, "\nnew chunk \n");
        printout(0, " idx  occ   sIx   pIx   cIx\n");
        for (int i = 0; i < this->maxNodesPerChunk; i++) {
            int status = this->nodeStackStatus[iShift + i];
            int sIdx = this->nodeChunks[iChunk][i].serialIx;
            int pIdx = this->nodeChunks[iChunk][i].parentSerialIx;
            int cIdx = this->nodeChunks[iChunk][i].childSerialIx;
            printout(0, std::setw(4) << n++);
            printout(0, std::setw(4) << status);
            printout(0, std::setw(6) << sIdx);
            printout(0, std::setw(6) << pIdx);
            printout(0, std::setw(6) << cIdx << "   ");
            if (status == 1) printout(0, this->nodeChunks[iChunk][i].getNodeIndex());
            printout(0, "\n");
        }
    }
}

template class FunctionNodeAllocator<1>;
template class FunctionNodeAllocator<2>;
template class FunctionNodeAllocator<3>;

} // namespace mrcpp
