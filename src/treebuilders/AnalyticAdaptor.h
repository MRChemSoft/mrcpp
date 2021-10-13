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

#pragma once

#include "MRCPP/constants.h"
#include "TreeAdaptor.h"

namespace mrcpp {

template <int D> class AnalyticAdaptor final : public TreeAdaptor<D> {
public:
    AnalyticAdaptor(const RepresentableFunction<D> &f, int ms)
            : TreeAdaptor<D>(ms)
            , func(&f) {}

private:
    const RepresentableFunction<D> *func;

    bool splitNode(const MWNode<D> &node) const override {
        int scale = node.getScale();
        int nQuadPts = node.getKp1();
        if (this->func->isVisibleAtScale(scale, nQuadPts)) return false;
        auto ub = node.getUpperBounds();
        auto lb = node.getLowerBounds();
        if (this->func->isZeroOnInterval(lb.data(), ub.data())) return false;
        return true;
    }
};

} // namespace mrcpp
