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

#include "OperatorNodeAllocator.h"
#include "OperatorNode.h"
#include "OperatorTree.h"
#include "utils/Printer.h"

namespace mrcpp {

/** SerialTree class constructor.
 * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
 * Initializes rootNodes to represent the zero function and allocate their nodes.
 * NOTES:
 * Serial trees are made of projected nodes, and include gennodes and loose nodes separately.
 * All created (using class creator) Projected nodes or GenNodes are loose nodes.
 * Loose nodes have their coeff in serial Tree, but not the node part.
 * Projected nodes and GenNodes that are created by their creator, are detroyed by destructor ~ProjectedNode and
 * ~GenNode. Serial tree nodes are not using the destructors, but explicitely call to deallocNodes or deallocGenNodes
 * Gen nodes and loose nodes are not counted with MWTree->[in/de]crementNodeCount()
 */
OperatorNodeAllocator::OperatorNodeAllocator(OperatorTree *tree)
        : NodeAllocator<2>(tree, nullptr)
        , lastNode(nullptr) {
    // reserve space for chunk pointers to avoid excessive reallocation
    this->nodeChunks.reserve(100);
    this->nodeCoeffChunks.reserve(100);

    this->topStack = 0;
    this->coeffsPerNode = 4 * this->tree_p->getKp1_d();
    this->maxNodesPerChunk = 1024;

    MRCPP_INIT_OMP_LOCK();
}

/** SerialTree destructor. */
OperatorNodeAllocator::~OperatorNodeAllocator() {
    for (auto &chunk : this->nodeChunks) delete[](char *) chunk;
    for (auto &chunk : this->nodeCoeffChunks) delete[] chunk;
    this->nodeStackStatus.clear();
    MRCPP_DESTROY_OMP_LOCK();
}

int OperatorNodeAllocator::getNChunks() const {
    return this->nodeChunks.size();
}

int OperatorNodeAllocator::getNodeChunkSize() const {
    return this->maxNodesPerChunk * sizeof(OperatorNode);
}

int OperatorNodeAllocator::getCoeffChunkSize() const {
    return this->maxNodesPerChunk * this->coeffsPerNode * sizeof(double);
}

OperatorNode * OperatorNodeAllocator::getNode_p(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int cIdx = sIdx % this->maxNodesPerChunk;   // position in chunk
    OperatorNode *node = this->nodeChunks[chunk] + cIdx;
    MRCPP_UNSET_OMP_LOCK();
    return node;
}

double * OperatorNodeAllocator::getCoef_p(int sIdx) {
    MRCPP_SET_OMP_LOCK();
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int idx = sIdx % this->maxNodesPerChunk;   // position in chunk
    double *coefs = this->nodeCoeffChunks[chunk] + idx * this->coeffsPerNode;
    MRCPP_UNSET_OMP_LOCK();
    return coefs;
}

OperatorNode * OperatorNodeAllocator::getNodeNoLock(int sIdx) {
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int cIdx = sIdx % this->maxNodesPerChunk;   // position in chunk
    OperatorNode *node = this->nodeChunks[chunk] + cIdx;
    return node;
}

double * OperatorNodeAllocator::getCoefNoLock(int sIdx) {
    int chunk = sIdx / this->maxNodesPerChunk; // which chunk
    int idx = sIdx % this->maxNodesPerChunk;   // position in chunk
    double *coefs = this->nodeCoeffChunks[chunk] + idx * this->coeffsPerNode;
    return coefs;
}

int OperatorNodeAllocator::alloc(int nAlloc, bool coeff) {
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

void OperatorNodeAllocator::dealloc(int serialIx) {
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

void OperatorNodeAllocator::appendChunk(bool coeff) {
    // make coeff chunk
    if (coeff) {
        double *c_chunk = new double[getCoeffChunkSize()];
        this->nodeCoeffChunks.push_back(c_chunk);
    }

    // make node chunk
    auto n_chunk = (OperatorNode *)new char[getNodeChunkSize()];
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

} // namespace mrcpp
