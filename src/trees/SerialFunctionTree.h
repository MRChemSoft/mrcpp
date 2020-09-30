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

#include <vector>

#include "SerialTree.h"
#include "utils/omp_utils.h"

namespace mrcpp {

template <int D> class SerialFunctionTree final : public SerialTree<D> {
public:
    SerialFunctionTree(FunctionTree<D> *tree, SharedMemory *sh_mem);
    SerialFunctionTree(const SerialFunctionTree<D> &tree) = delete;
    SerialFunctionTree<D> &operator=(const SerialFunctionTree<D> &tree) = delete;
    ~SerialFunctionTree() override;

    void allocRoots(MWTree<D> &tree) override;
    void allocChildren(MWNode<D> &parent) override;
    void allocGenChildren(MWNode<D> &parent) override;

    void deallocNodes(int serialIx) override;
    void deallocGenNodes(int serialIx) override;
    void deallocGenNodeChunks() override;

    int getNChunks() const { return this->nodeChunks.size(); }
    int getNChunksUsed() const;

    int shrinkChunks();

    std::vector<ProjectedNode<D> *> nodeChunks;
    std::vector<double *> nodeCoeffChunks;

    ProjectedNode<D> *sNodes; // serial ProjectedNodes
    GenNode<D> *sGenNodes;    // serial GenNodes

    std::vector<GenNode<D> *> genNodeChunks;
    std::vector<double *> genNodeCoeffChunks;

    int nGenNodes; // number of GenNodes already defined

    double **genCoeffStack;

    //    int *genNodeStackStatus;
    std::vector<int> genNodeStackStatus;

    void rewritePointers();

    void clear(int n);

protected:
    int maxGenNodes;      // max number of Gen nodes that can be defined
    int sizeGenNodeCoeff; // size of coeff for one Gen node

    char *cvptr_ProjectedNode; // virtual table pointer for ProjectedNode
    char *cvptr_GenNode;       // virtual table pointer for GenNode

    ProjectedNode<D> *lastNode; // pointer just after the last active node, i.e. where to put next node
    GenNode<D> *lastGenNode;    // pointer just after the last active Gen node, i.e. where to put next node

    ProjectedNode<D> *allocNodes(int nAlloc, int *serialIx, double **coefs_p);
    GenNode<D> *allocGenNodes(int nAlloc, int *serialIx, double **coefs_p);

private:
#ifdef MRCPP_HAS_OMP
    omp_lock_t omp_lock;
#endif
};

} // namespace mrcpp
