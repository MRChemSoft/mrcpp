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

#include "NodeAllocator.h"

#include <stack>

#include "MWNode.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "OperatorTree.h"
#include "OperatorNode.h"

#include "utils/Printer.h"
#include "utils/mpi_utils.h"

namespace mrcpp {

template <int D> NodeAllocator<D>::NodeAllocator(FunctionTree<D> *tree, SharedMemory *mem, bool gen)
        : tree_p(tree)
        , shmem_p(mem) {
    // reserve space for chunk pointers to avoid excessive reallocation
    this->nodeChunks.reserve(100);
    this->coeffChunks.reserve(100);

    int tDim = (gen) ? 1 : (1 << D); // genNodes have only one block
    this->coeffsPerNode = tDim * getTree().getKp1_d();

    if (D < 3) {
        // define from number of nodes per chunk
        this->maxNodesPerChunk = 64;
    } else {
        // 2 MB small for no waisting place, but large enough so that latency and overhead work is negligible
        int sizePerChunk = 2 * 1024 * 1024;
        this->maxNodesPerChunk = (sizePerChunk / this->coeffsPerNode / sizeof(double) / 8) * 8;
    }

    FunctionNode<D> tmp;
    this->cvptr = *(char **)(&tmp);
    this->sizeOfNode = sizeof(FunctionNode<D>);

    MRCPP_INIT_OMP_LOCK();
}

template <> NodeAllocator<2>::NodeAllocator(OperatorTree *tree, SharedMemory *mem)
        : tree_p(tree)
        , shmem_p(mem) {
    // reserve space for chunk pointers to avoid excessive reallocation
    this->nodeChunks.reserve(100);
    this->coeffChunks.reserve(100);

    this->coeffsPerNode = 4 * getTree().getKp1_d();
    this->maxNodesPerChunk = 1024;

    OperatorNode tmp;
    this->cvptr = *(char **)(&tmp);
    this->sizeOfNode = sizeof(OperatorNode);

    MRCPP_INIT_OMP_LOCK();
}

template <int D> NodeAllocator<D>::NodeAllocator(OperatorTree *tree, SharedMemory *mem) {
    NOT_REACHED_ABORT;
}

template <int D> NodeAllocator<D>::~NodeAllocator() {
    for (auto &chunk : this->nodeChunks) delete[](char *) chunk;
    if (not isShared()) // if the data is shared, it must be freed by MPI_Win_free
        for (auto &chunk : this->coeffChunks) delete[] chunk;
    this->stackStatus.clear();
    MRCPP_DESTROY_OMP_LOCK();
}

/** reset the start node counter */
template <int D> void NodeAllocator<D>::clear(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    std::fill(this->stackStatus.begin() + sIdx, this->stackStatus.end(), 0);
    this->topStack = sIdx;
    this->last_p = getNodeNoLock(sIdx);

    if (isShared())
        getMemory().sh_end_ptr = getMemory().sh_start_ptr + (sIdx / this->maxNodesPerChunk + 1) * this->coeffsPerNode * this->maxNodesPerChunk;

    MRCPP_UNSET_OMP_LOCK();
}

template <int D> MWNode<D> * NodeAllocator<D>::getNode_p(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int cIdx = sIdx % this->maxNodesPerChunk;   // position in chunk
    MWNode<D> *node = this->nodeChunks[chunk] + cIdx;
    MRCPP_UNSET_OMP_LOCK();
    return node;
}

template <int D> double * NodeAllocator<D>::getCoef_p(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int idx = sIdx % this->maxNodesPerChunk;   // position in chunk
    double *coefs = this->coeffChunks[chunk] + idx * this->coeffsPerNode;
    MRCPP_UNSET_OMP_LOCK();
    return coefs;
}

template <int D> MWNode<D> * NodeAllocator<D>::getNodeNoLock(int sIdx) {
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int cIdx = sIdx % this->maxNodesPerChunk;   // position in chunk
    MWNode<D> *node = this->nodeChunks[chunk] + cIdx;
    return node;
}

template <int D> double * NodeAllocator<D>::getCoefNoLock(int sIdx) {
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int idx = sIdx % this->maxNodesPerChunk;   // position in chunk
    double *coefs = this->coeffChunks[chunk] + idx * this->coeffsPerNode;
    return coefs;
}

template <int D> int NodeAllocator<D>::alloc(int nAlloc, bool coeff) {
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
    auto &status = this->stackStatus;
    for (int i = sIdx; i < sIdx + nAlloc; i++) {
        if (status[i] != 0) MSG_ERROR(" NodeStackStatus: not available [" << i << "] : " << status[i]);
        status[i] = 1;
    }

    // advance stack pointers
    this->nNodes += nAlloc;
    this->topStack += nAlloc;
    this->last_p = getNodeNoLock(sIdx) + nAlloc;
    MRCPP_UNSET_OMP_LOCK();

    return sIdx;
}

template <int D> void NodeAllocator<D>::dealloc(int serialIx) {
    MRCPP_SET_OMP_LOCK();
    this->stackStatus[serialIx] = 0; // mark as available
    if (serialIx == this->topStack - 1) {  // top of stack
        while (this->stackStatus[this->topStack - 1] == 0) {
            this->topStack--;
            if (this->topStack < 1) break;
        }
        // has to redefine last_p
        this->last_p = getNodeNoLock(this->topStack);
    }
    this->nNodes--;
    MRCPP_UNSET_OMP_LOCK();
}

template <int D> void NodeAllocator<D>::initChunk(int iChunk, bool coeff) {
    MRCPP_SET_OMP_LOCK();
    if (iChunk >= getNChunks()) appendChunk(coeff);
    MRCPP_UNSET_OMP_LOCK();
}

template <int D> void NodeAllocator<D>::appendChunk(bool coeff) {
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
        this->coeffChunks.push_back(c_chunk);
    }

    // make node chunk
    auto n_chunk = (MWNode<D> *)new char[getNodeChunkSize()];
    for (int i = 0; i < this->maxNodesPerChunk; i++) {
        n_chunk[i].serialIx = -1;
        n_chunk[i].parentSerialIx = -1;
        n_chunk[i].childSerialIx = -1;
    }
    this->nodeChunks.push_back(n_chunk);

    // append to stackStatus
    int oldsize = this->stackStatus.size();
    int newsize = oldsize + this->maxNodesPerChunk;
    this->stackStatus.resize(newsize);
    std::fill(this->stackStatus.begin() + oldsize, this->stackStatus.end(), 0);
}

/** Fill all holes in the chunks with occupied nodes, then remove all empty chunks */
template <int D> int NodeAllocator<D>::compress() {
    MRCPP_SET_OMP_LOCK();
    int nAlloc = (1 << D);
    if (this->maxNodesPerChunk * this->nodeChunks.size() <=
        getTree().getNNodes() + this->maxNodesPerChunk + nAlloc - 1) {
        MRCPP_UNSET_OMP_LOCK();
        return 0; // no chunks to remove
    }

    int posocc = 0;
    int posavail = getTree().getRootBox().size(); // start after root nodes
    while (true) {
        posavail = findNextAvailable(posavail, nAlloc);
        if (posavail >= this->topStack) break; // treated all nodes

        posocc = findNextOccupied(posavail);
        if (posocc >= this->topStack) break; // treated all nodes

        moveNodes(nAlloc, posocc, posavail);
    }

    // find the last used node
    posocc = this->topStack - 1;
    while (this->stackStatus[posocc] == 0 and posocc > 0) posocc--;
    this->topStack = posocc + 1;
    this->last_p = getNodeNoLock(this->topStack);

    int nChunksDeleted = deleteUnusedChunks();
    getTree().resetEndNodeTable();

    MRCPP_UNSET_OMP_LOCK();
    return nChunksDeleted;
}

template <int D> int NodeAllocator<D>::deleteUnusedChunks() {
    // number of occupied chunks
    int nChunksTotal = getNChunks();
    int nChunksUsed = getNChunksUsed();
    for (int i = nChunksUsed; i < nChunksTotal; i++) delete[](char *)(this->nodeChunks[i]);

    if (isShared()) {
        // shared coefficients cannot be fully deallocated, only pointer is moved.
       getMemory().sh_end_ptr -= (nChunksTotal - nChunksUsed) * this->coeffsPerNode * this->maxNodesPerChunk;
    } else {
        for (int i = nChunksUsed; i < nChunksTotal; i++) delete[] this->coeffChunks[i];
    }

    // shrink the stacks
    this->nodeChunks.resize(nChunksUsed);
    this->coeffChunks.resize(nChunksUsed);
    this->stackStatus.resize(nChunksUsed * this->maxNodesPerChunk);
    return nChunksTotal - nChunksUsed;
}

template <int D> void NodeAllocator<D>::moveNodes(int nNodes, int srcIdx, int dstIdx) {
    auto *srcNode = getNodeNoLock(srcIdx);
    auto *dstNode = getNodeNoLock(dstIdx);

    // check that all siblings are consecutive. Should never be root node.
    for (int i = 0; i < nNodes; i++) assert(this->stackStatus[dstIdx + i] == 0);
    for (int i = 1; i < nNodes; i++) assert((srcNode + i)->parent->serialIx == srcNode->parent->serialIx); // siblings

    // just copy everything "as is"
    for (int i = 0; i < nNodes * this->sizeOfNode; i++) ((char *)dstNode)[i] = ((char *)srcNode)[i];

    // coefs have new adresses
    double *coefs_p = getCoefNoLock(dstIdx);
    for (int i = 0; i < nNodes; i++) (dstNode + i)->coefs = coefs_p + i * getNCoefs();

    // copy coefs to new adress
    if (not isShared()) {
        for (int i = 0; i < nNodes * this->coeffsPerNode; i++) dstNode->coefs[i] = srcNode->coefs[i];
    } else {
        if (getMemory().rank == 0) // only master copy the data. careful with sync
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
    for (int i = 0; i < nNodes; i++) this->stackStatus[dstIdx + i] = 1;
    dstIdx += nNodes;

    // delete "old" nodes
    for (int i = 0; i < nNodes; i++) this->stackStatus[srcIdx + i] = 0;
    for (int i = 0; i < nNodes; i++) (srcNode + i)->serialIx = -1;
}

template <int D> int NodeAllocator<D>::findNextAvailable(int pos, int nAlloc) const {
    // Last positions on a chunk cannot be used if there is no place for nAlloc siblings on the same chunk
    bool chunkTooSmall = (pos + nAlloc - 1) / this->maxNodesPerChunk != pos / this->maxNodesPerChunk;
    bool posIsOccupied = (this->stackStatus[pos] != 0);
    bool endOfStack = (pos >= this->topStack);
    while ((posIsOccupied or chunkTooSmall) and not(endOfStack)) {
        pos++;
        chunkTooSmall = (pos + nAlloc - 1) / this->maxNodesPerChunk != pos / this->maxNodesPerChunk;
        posIsOccupied = (this->stackStatus[pos] != 0);
        endOfStack = (pos >= this->topStack);
    }
    return pos;
}

template <int D> int NodeAllocator<D>::findNextOccupied(int pos) const {
    bool posIsAvailable = (this->stackStatus[pos] == 0);
    bool endOfStack = (pos >= this->topStack);
    while (posIsAvailable and not(endOfStack)) {
        pos++;
        posIsAvailable = (this->stackStatus[pos] == 0);
        endOfStack = (pos >= this->topStack);
    }
    return pos;
}

/** Traverse tree and redefine pointer, counter and tables. */
template <int D> void NodeAllocator<D>::rewritePointers(bool coeff) {
    MRCPP_SET_OMP_LOCK();
    this->nNodes = 0;
    getTree().nodesAtDepth.clear();
    getTree().squareNorm = 0.0;
    getTree().clearEndNodeTable();

    // reinitialize stacks
    // nodeChunks have been adapted to receiving tree and maybe larger than stackStatus
    int nodeCount = this->nodeChunks.size() * this->maxNodesPerChunk;
    this->stackStatus.resize(nodeCount);
    std::fill(this->stackStatus.begin(), this->stackStatus.end(), 0);

    NodeBox<D> &rootbox = getTree().getRootBox();
    MWNode<D> **roots = rootbox.getNodes();

    std::stack<MWNode<D> *> stack;
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
        getTree().incrementNodeCount(node->getScale());
        if (node->isEndNode()) getTree().squareNorm += node->getSquareNorm();
        if (node->isEndNode()) getTree().endNodeTable.push_back(node);

        // normally (intel) the virtual table does not change, but we overwrite anyway
        *(char **)(node) = this->cvptr;

        node->tree = this->tree_p;
        node->coefs = (coeff) ? getCoefNoLock(sIdx) : nullptr;
        node->parent = (pIdx >= 0) ? getNodeNoLock(pIdx) : nullptr;

        stack.pop();
        auto *child_p = getNodeNoLock(cIdx);
        for (int i = 0; i < node->getNChildren(); i++) {
            node->children[i] = child_p;
            stack.push(child_p);
            child_p++;
        }
        this->stackStatus[sIdx] = 1; // occupied
    }
    this->last_p = getNodeNoLock(this->topStack);
    MRCPP_UNSET_OMP_LOCK();
}

template <int D> void NodeAllocator<D>::print() const {
    int n = 0;
    for (int iChunk = 0; iChunk < getNChunks(); iChunk++) {
        int iShift = iChunk * this->maxNodesPerChunk;
        printout(0, "\nnew chunk \n");
        printout(0, " idx  occ   sIx   pIx   cIx\n");
        for (int i = 0; i < this->maxNodesPerChunk; i++) {
            int status = this->stackStatus[iShift + i];
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

template class NodeAllocator<1>;
template class NodeAllocator<2>;
template class NodeAllocator<3>;

} // namespace mrcpp
