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

#pragma once

#include "TreeAdaptor.h"
#include "trees/FunctionTreeVector.h"
#include "utils/Printer.h"

namespace mrcpp {

template <int D> class MultiplicationAdaptor : public TreeAdaptor<D> {
public:
    MultiplicationAdaptor(double pr, int ms, FunctionTreeVector<D> &t)
            : TreeAdaptor<D>(ms)
            , prec(pr)
            , trees(t) {}
    ~MultiplicationAdaptor() override = default;

protected:
    double prec;
    mutable FunctionTreeVector<D> trees;

    bool splitNode(const MWNode<D> &node) const override {
        if (this->trees.size() != 2) MSG_ERROR("Invalid tree vec size: " << this->trees.size());
        auto multPrec = 1.0;
        auto &pNode0 = get_func(trees, 0).getNode(node.getNodeIndex());
        auto &pNode1 = get_func(trees, 1).getNode(node.getNodeIndex());
        double maxW0 = std::sqrt(pNode0.getMaxWSquareNorm());
        double maxW1 = std::sqrt(pNode1.getMaxWSquareNorm());
        double maxS0 = std::sqrt(pNode0.getMaxSquareNorm());
        double maxS1 = std::sqrt(pNode1.getMaxSquareNorm());
        if (pNode0.isGenNode()) {
            maxW0 = 0.0;
            maxS0 = std::sqrt(std::pow(2.0, D * node.getScale()) * pNode0.getSquareNorm());
        }
        if (pNode1.isGenNode()) {
            maxW1 = 0.0;
            maxS1 = std::sqrt(std::pow(2.0, D * node.getScale()) * pNode1.getSquareNorm());
        }
        // The wavelet contribution (in the product of node0 and node1) can be approximated as
        multPrec = maxW0 * maxS1 + maxW1 * maxS0 + maxW0 * maxW1;

        // Note: this never refine deeper than one scale more than input tree grids, because when wavelets are zero
        // for both input trees, multPrec=0 In addition, we force not to refine deeper than input tree grids
        if (multPrec > this->prec and not(pNode0.isLeafNode() and pNode1.isLeafNode())) {
            return true;
        } else {
            return false;
        }
    }
};

} // namespace mrcpp
