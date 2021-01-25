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

#include <Eigen/Core>
#include <vector>

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

template <int D> class SerialTree {
public:
    SerialTree(MWTree<D> *tree, SharedMemory *mem);
    SerialTree(const SerialTree<D> &tree) = delete;
    SerialTree<D> &operator=(const SerialTree<D> &tree) = delete;
    virtual ~SerialTree() = default;

    MWTree<D> *getTree() { return this->tree_p; }
    SharedMemory *getMemory() { return this->shMem; }

    bool isShared() const {
        if (this->shMem == nullptr) return false;
        return true;
    }

    virtual void allocRoots(MWTree<D> &tree) = 0;
    virtual void allocChildren(MWNode<D> &parent) = 0;
    virtual void allocChildrenNoCoeff(MWNode<D> &parent) = 0;
    virtual void allocGenChildren(MWNode<D> &parent) = 0;

    virtual void deallocNodes(int serialIx) = 0;
    virtual void deallocGenNodes(int serialIx) = 0;
    virtual void deallocGenNodeChunks() = 0;

    void S_mwTransform(double *coeff_in, double *coeff_out, bool readOnlyScaling, int stride, bool overwrite = true);
    void S_mwTransformBack(double *coeff_in, double *coeff_out, int stride);

    int nNodes; // number of Nodes already defined
    int maxNodesPerChunk;
    std::vector<int> nodeStackStatus;
    int sizeNodeCoeff; // size of coeff for one node
    double **coeffStack;
    int maxNodes; // max number of nodes that can be defined

protected:
    MWTree<D> *tree_p;
    SharedMemory *shMem;
};

} // namespace mrcpp
