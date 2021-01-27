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

template <int D> class NodeAllocator {
public:
    NodeAllocator(MWTree<D> *tree, SharedMemory *mem) : tree_p(tree), shMem(mem) {}
    NodeAllocator(const NodeAllocator<D> &tree) = delete;
    NodeAllocator<D> &operator=(const NodeAllocator<D> &tree) = delete;
    virtual ~NodeAllocator() = default;

    MWTree<D> *getTree() { return this->tree_p; }
    SharedMemory *getMemory() { return this->shMem; }

    bool isShared() const { return (this->shMem != nullptr); }

    virtual void allocRoots(MWTree<D> &tree) = 0;
    virtual void allocChildren(MWNode<D> &parent) = 0;
    virtual void allocChildrenNoCoeff(MWNode<D> &parent) = 0;
    virtual void deallocNodes(int serialIx) = 0;

    virtual int getNChunks() const = 0;

protected:
    int nNodes{0};                  // number of Nodes already defined
    int maxNodes{0};                // max number of nodes that can be defined
    int coeffsPerNode{0};           // number of coeff for one node
    int maxNodesPerChunk{0};        // max number of nodes per allocation
    std::vector<int> nodeStackStatus;
    std::vector<double *> nodeCoeffChunks;

    MWTree<D> *tree_p{nullptr};     // pointer to external object
    SharedMemory *shMem{nullptr}; // pointer to external object

#ifdef MRCPP_HAS_OMP
    omp_lock_t omp_lock;
#endif
};


} // namespace mrcpp
