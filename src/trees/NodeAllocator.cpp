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
    if (sIdx < 0 or sIdx >= this->stackStatus.size()) MSG_ABORT("Invalid serial index: " << sIdx);
    std::fill(this->stackStatus.begin() + sIdx, this->stackStatus.end(), 0);
    this->topStack = sIdx;
    this->last_p = getNodeNoLock(sIdx);

    if (isShared())
        getMemory().sh_end_ptr = getMemory().sh_start_ptr + (sIdx / this->maxNodesPerChunk + 1) * this->coeffsPerNode * this->maxNodesPerChunk;

    MRCPP_UNSET_OMP_LOCK();
}

template <int D> MWNode<D> * NodeAllocator<D>::getNode_p(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    auto *node = getNodeNoLock(sIdx);
    MRCPP_UNSET_OMP_LOCK();
    return node;
}

template <int D> double * NodeAllocator<D>::getCoef_p(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    auto *coefs = getCoefNoLock(sIdx);
    MRCPP_UNSET_OMP_LOCK();
    return coefs;
}

template <int D> MWNode<D> * NodeAllocator<D>::getNodeNoLock(int sIdx) {
    if (sIdx < 0 or sIdx >= this->stackStatus.size()) return nullptr;
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int cIdx = sIdx % this->maxNodesPerChunk;  // position in chunk
    return this->nodeChunks[chunk] + cIdx;
}

template <int D> double * NodeAllocator<D>::getCoefNoLock(int sIdx) {
    if (sIdx < 0 or sIdx >= this->stackStatus.size()) return nullptr;
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int idx = sIdx % this->maxNodesPerChunk;   // position in chunk
    return this->coeffChunks[chunk] + idx * this->coeffsPerNode;
}

template <int D> int NodeAllocator<D>::alloc(int nNodes, bool coeff) {
    MRCPP_SET_OMP_LOCK();
    if (nNodes <= 0 or nNodes > this->maxNodesPerChunk) MSG_ABORT("Cannot allocate " << nNodes << " nodes");

    // move topstack to start of next chunk if current chunk is too small
    int cIdx = this->topStack % (this->maxNodesPerChunk);
    bool chunkOverflow = ((cIdx + nNodes) > this->maxNodesPerChunk);
    if (chunkOverflow) this->topStack = this->maxNodesPerChunk * ((this->topStack + nNodes - 1) / this->maxNodesPerChunk);

    // append chunk if necessary
    int chunk = this->topStack / this->maxNodesPerChunk;
    bool needNewChunk = (chunk >= this->nodeChunks.size());
    if (needNewChunk) appendChunk(coeff);

    // return value is index of first new node
    auto sIdx = this->topStack;

    // fill stack status
    auto &status = this->stackStatus;
    for (int i = sIdx; i < sIdx + nNodes; i++) {
        if (status[i] != 0) MSG_ERROR(" NodeStackStatus: not available [" << i << "] : " << status[i]);
        status[i] = 1;
    }

    // advance stack pointers
    this->nNodes += nNodes;
    this->topStack += nNodes;
    this->last_p = getNodeNoLock(sIdx) + nNodes;
    MRCPP_UNSET_OMP_LOCK();

    return sIdx;
}

template <int D> void NodeAllocator<D>::dealloc(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    if (sIdx < 0 or sIdx >= this->stackStatus.size()) MSG_ABORT("Invalid serial index: " << sIdx);
    this->stackStatus[sIdx] = 0; // mark as available
    if (sIdx == this->topStack - 1) {  // top of stack
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

template <int D> void NodeAllocator<D>::init(int nChunks, bool coefs) {
    MRCPP_SET_OMP_LOCK();
    if (nChunks <= 0) MSG_ABORT("Invalid number of chunks: " << nChunks);
    for (int i = getNChunks(); i < nChunks; i++) appendChunk(coefs);

    // reinitialize stacks
    int nodeCount = this->nodeChunks.size() * this->maxNodesPerChunk;
    this->stackStatus.resize(nodeCount);
    std::fill(this->stackStatus.begin(), this->stackStatus.end(), 0);
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
    int nNodes = (1 << D);
    if (this->maxNodesPerChunk * this->nodeChunks.size() <=
        getTree().getNNodes() + this->maxNodesPerChunk + nNodes - 1) {
        MRCPP_UNSET_OMP_LOCK();
        return 0; // nothing to compress
    }

    int posocc = 0;
    int posavail = getTree().getRootBox().size(); // start after root nodes
    while (true) {
        posavail = findNextAvailable(posavail, nNodes);
        if (posavail >= this->topStack) break; // treated all nodes

        posocc = findNextOccupied(posavail);
        if (posocc >= this->topStack) break; // treated all nodes

        moveNodes(nNodes, posocc, posavail);
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
    assert(nChunksTotal >= nChunksUsed);
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
    assert(nNodes > 0);
    assert(nNodes <= this->maxNodesPerChunk);

    auto *srcNode = getNodeNoLock(srcIdx);
    auto *dstNode = getNodeNoLock(dstIdx);
    assert(srcNode != nullptr);
    assert(dstNode != nullptr);

    // check that all siblings are consecutive. Should never be root node.
    for (int i = 0; i < nNodes; i++) assert(this->stackStatus[dstIdx + i] == 0);
    for (int i = 1; i < nNodes; i++) assert((srcNode + i)->parent->serialIx == srcNode->parent->serialIx); // siblings

    // just copy everything "as is"
    for (int i = 0; i < nNodes * this->sizeOfNode; i++) ((char *)dstNode)[i] = ((char *)srcNode)[i];

    // coefs have new adresses
    double *coefs_p = getCoefNoLock(dstIdx);
    if (coefs_p == nullptr) NOT_IMPLEMENTED_ABORT; // Nodes without coefs not handled atm
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

// Last positions on a chunk cannot be used if there is no place for nNodes siblings on the same chunk
template <int D> int NodeAllocator<D>::findNextAvailable(int sIdx, int nNodes) const {
    assert(sIdx >= 0);
    assert(sIdx < this->stackStatus.size());
    assert(nNodes >= 0);
    assert(nNodes < this->maxNodesPerChunk);
    bool chunkTooSmall = (sIdx + nNodes - 1) / this->maxNodesPerChunk != sIdx / this->maxNodesPerChunk;
    bool posIsOccupied = (this->stackStatus[sIdx] != 0);
    bool endOfStack = (sIdx >= this->topStack);
    while ((posIsOccupied or chunkTooSmall) and not(endOfStack)) {
        sIdx++;
        chunkTooSmall = (sIdx + nNodes - 1) / this->maxNodesPerChunk != sIdx / this->maxNodesPerChunk;
        posIsOccupied = (this->stackStatus[sIdx] != 0);
        endOfStack = (sIdx >= this->topStack);
    }
    assert(sIdx >= 0);
    assert(sIdx < this->stackStatus.size());
    return sIdx;
}

template <int D> int NodeAllocator<D>::findNextOccupied(int sIdx) const {
    assert(sIdx >= 0);
    assert(sIdx < this->stackStatus.size());
    bool posIsAvailable = (this->stackStatus[sIdx] == 0);
    bool endOfStack = (sIdx >= this->topStack);
    while (posIsAvailable and not(endOfStack)) {
        sIdx++;
        posIsAvailable = (this->stackStatus[sIdx] == 0);
        endOfStack = (sIdx >= this->topStack);
    }
    assert(sIdx >= 0);
    assert(sIdx < this->stackStatus.size());
    return sIdx;
}

/** Traverse tree and redefine pointer, counter and tables. */
template <int D> void NodeAllocator<D>::reassemble() {
    MRCPP_SET_OMP_LOCK();
    this->nNodes = 0;
    getTree().nodesAtDepth.clear();
    getTree().squareNorm = 0.0;
    getTree().clearEndNodeTable();

    NodeBox<D> &rootbox = getTree().getRootBox();
    MWNode<D> **roots = rootbox.getNodes();

    std::stack<MWNode<D> *> stack;
    for (int rIdx = 0; rIdx < rootbox.size(); rIdx++) {
        auto *root_p = getNodeNoLock(rIdx);
        assert(root_p != nullptr);
        stack.push(root_p);
        roots[rIdx] = root_p;
    }
    this->topStack = 0;
    while (not stack.empty()) {
        auto *node_p = stack.top();
        assert(node_p != nullptr);
        auto sIdx = node_p->serialIx;
        auto pIdx = node_p->parentSerialIx;
        auto cIdx = node_p->childSerialIx;

        this->nNodes++;
        this->topStack = std::max(this->topStack, sIdx + 1);
        getTree().incrementNodeCount(node_p->getScale());
        if (node_p->isEndNode()) getTree().squareNorm += node_p->getSquareNorm();
        if (node_p->isEndNode()) getTree().endNodeTable.push_back(node_p);

        // normally (intel) the virtual table does not change, but we overwrite anyway
        *(char **)(node_p) = this->cvptr;

        node_p->tree = this->tree_p;
        node_p->coefs = getCoefNoLock(sIdx);
        node_p->parent = getNodeNoLock(pIdx);

        stack.pop();
        auto *child_p = getNodeNoLock(cIdx);
        for (int i = 0; i < node_p->getNChildren(); i++) {
            assert(child_p != nullptr);
            node_p->children[i] = child_p;
            stack.push(child_p);
            child_p++;
        }
        this->stackStatus[sIdx] = 1; // occupied
    }
    this->last_p = getNodeNoLock(this->topStack);
    assert(this->last_p != nullptr);
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
