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

/**
 *
 *
 *  \date Jul, 2016
 *  \author Peter Wind <peter.wind@uit.no> \n
 *  CTCC, University of Tromsø
 *
 */

#pragma once

#include "NodeAllocator.h"
#include "FunctionNode.h"

namespace mrcpp {

template <int D> class FunctionNodeAllocator final : public NodeAllocator<D> {
public:
    FunctionNodeAllocator(FunctionTree<D> *tree, SharedMemory *sh_mem, bool gen = false);
    FunctionNodeAllocator(const FunctionNodeAllocator<D> &tree) = delete;
    FunctionNodeAllocator<D> &operator=(const FunctionNodeAllocator<D> &tree) = delete;
    ~FunctionNodeAllocator() override;

    int getNChunks() const override { return this->nodeChunks.size(); }
    int getNChunksUsed() const { return (this->topStack + this->maxNodesPerChunk - 1) / this->maxNodesPerChunk; }
    int getNodeChunkSize() const { return this->maxNodesPerChunk * sizeof(FunctionNode<D>); }
    int getCoeffChunkSize() const { return this->maxNodesPerChunk * this->coeffsPerNode * sizeof(double); }

    double *getCoeffChunk(int i) { return this->nodeCoeffChunks[i]; }
    FunctionNode<D> *getNodeChunk(int i) { return this->nodeChunks[i]; }

    int alloc(int nAlloc, bool coeff = true);
    void dealloc(int serialIx) override;

    void initChunk(int iChunk, bool coeff = true);
    int shrinkChunks();

    void rewritePointers(bool coeff = true);
    void clear(int n);
    void print() const;

    double * getCoef_p(int sIdx);
    FunctionNode<D> * getNode_p(int sIdx);

protected:
    bool genNode{false};
    FunctionNode<D> *lastNode{nullptr}; // pointer just after the last active node, i.e. where to put next node
    std::vector<FunctionNode<D> *> nodeChunks{};

    double * getCoefNoLock(int sIdx);
    FunctionNode<D> * getNodeNoLock(int sIdx);
    void appendChunk(bool coeff = true);
};

} // namespace mrcpp
