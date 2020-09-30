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

class SerialOperatorTree final : public SerialTree<2> {
public:
    SerialOperatorTree(OperatorTree *tree);
    SerialOperatorTree(const SerialOperatorTree &tree) = delete;
    SerialOperatorTree &operator=(const SerialOperatorTree &tree) = delete;
    ~SerialOperatorTree() override;

    void allocRoots(MWTree<2> &tree) override;
    void allocChildren(MWNode<2> &parent) override;
    void allocGenChildren(MWNode<2> &parent) override;

    void deallocNodes(int serialIx) override;
    void deallocGenNodes(int serialIx) override;
    void deallocGenNodeChunks() override;

protected:
    OperatorNode *sNodes; // serial OperatorNodes

    std::vector<OperatorNode *> nodeChunks;
    std::vector<double *> nodeCoeffChunks;

    char *cvptr_OperatorNode; // virtual table pointer for OperatorNode
    OperatorNode *lastNode;   // pointer to the last active node

    OperatorNode *allocNodes(int nAlloc, int *serialIx, double **coefs_p);

private:
#ifdef MRCPP_HAS_OMP
    omp_lock_t omp_lock;
#endif
};

} // namespace mrcpp
