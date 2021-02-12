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
#include "utils/Printer.h"
#include "utils/tree_utils.h"

namespace mrcpp {

template <int D> class WaveletAdaptor : public TreeAdaptor<D> {
public:
    WaveletAdaptor(double pr, int ms, bool ap = false, double sf = 1.0)
            : TreeAdaptor<D>(ms)
            , absPrec(ap)
            , prec(pr)
            , splitFac(sf) {}
    ~WaveletAdaptor() override = default;

    void setPrecFunction(const std::function<double(const NodeIndex<D> &idx)> &prec_func) {
        this->precFunc = prec_func;
    }

protected:
    bool absPrec;
    double prec;
    double splitFac;
    std::function<double(const NodeIndex<D> &idx)> precFunc = [](const NodeIndex<D> &idx) { return 1.0; };

    bool splitNode(const MWNode<D> &node) const override {
        auto precFac = this->precFunc(node.getNodeIndex()); // returns 1.0 by default
        return tree_utils::split_check(node, this->prec * precFac, this->splitFac, this->absPrec);
    }
};

} // namespace mrcpp
