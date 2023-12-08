/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

/**
 *
 *  \date July, 2016
 *  \author Peter Wind <peter.wind@uit.no> \n
 *  CTCC, University of Tromsø
 *
 */

#pragma once

#include <vector>

#include "MRCPP/mrcpp_declarations.h"
#include "utils/omp_utils.h"

namespace mrcpp {

template <int D> class NodeAllocator final {
public:
    NodeAllocator(OperatorTree *tree, SharedMemory *mem, int coefsPerNode, int nodesPerChunk);
    NodeAllocator(FunctionTree<D> *tree, SharedMemory *mem, int coefsPerNode, int nodesPerChunk);
    NodeAllocator(const NodeAllocator<D> &tree) = delete;
    NodeAllocator<D> &operator=(const NodeAllocator<D> &tree) = delete;
    ~NodeAllocator();

    int alloc(int nNodes, bool coefs = true);
    void dealloc(int sIdx);
    void deallocAllCoeff();

    void init(int nChunks, bool coefs = true);

    int compress();
    void reassemble();
    int deleteUnusedChunks();

    int getNNodes() const { return this->nNodes; }
    int getNCoefs() const { return this->coefsPerNode; }
    int getNChunks() const { return this->nodeChunks.size(); }
    int getNChunksUsed() const { return (this->topStack + this->maxNodesPerChunk - 1) / this->maxNodesPerChunk; }
    int getNodeChunkSize() const { return this->maxNodesPerChunk * this->sizeOfNode; }
    int getCoefChunkSize() const { return this->maxNodesPerChunk * this->coefsPerNode * sizeof(double); }

    double * getCoef_p(int sIdx);
    MWNode<D> * getNode_p(int sIdx);

    double * getCoefChunk(int i) { return this->coefChunks[i]; }
    MWNode<D> * getNodeChunk(int i) { return this->nodeChunks[i]; }

    void print() const;

protected:
    int nNodes{0};                  // number of nodes actually in use
    int topStack{0};                // index of last node on stack
    int sizeOfNode{0};              // sizeof(NodeType)
    int coefsPerNode{0};            // number of coef for one node
    int maxNodesPerChunk{0};        // max number of nodes per allocation

    std::vector<int> stackStatus{};
    std::vector<double *> coefChunks{};
    std::vector<MWNode<D> *> nodeChunks{};

    char *cvptr{nullptr};           // pointer to virtual table
    MWNode<D> *last_p{nullptr};     // pointer just after the last active node, i.e. where to put next node
    MWTree<D> *tree_p{nullptr};     // pointer to external object
    SharedMemory *shmem_p{nullptr}; // pointer to external object

    bool isShared() const { return (this->shmem_p != nullptr); }
    MWTree<D> &getTree() { return *this->tree_p; }
    SharedMemory &getMemory() { return *this->shmem_p; }

    double * getCoefNoLock(int sIdx);
    MWNode<D> * getNodeNoLock(int sIdx);

    void moveNodes(int nNodes, int srcIdx, int dstIdx);
    void appendChunk(bool coefs);

    int findNextAvailable(int sIdx, int nNodes) const;
    int findNextOccupied(int sIdx) const;

#ifdef MRCPP_HAS_OMP
    omp_lock_t omp_lock;
#endif
};


} // namespace mrcpp
