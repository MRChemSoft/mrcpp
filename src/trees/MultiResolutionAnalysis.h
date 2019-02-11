/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "BoundingBox.h"
#include "core/MWFilter.h"
#include "core/ScalingBasis.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

template <int D> class MultiResolutionAnalysis final {
public:
    MultiResolutionAnalysis(const BoundingBox<D> &bb, const ScalingBasis &sb, int depth = MaxDepth);
    MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra);
    MultiResolutionAnalysis &operator=(const MultiResolutionAnalysis &mra) = delete;

    int getOrder() const { return this->basis.getScalingOrder(); }
    int getMaxDepth() const { return this->maxDepth; }
    int getMaxScale() const { return this->world.getScale() + this->maxDepth; }

    const MWFilter &getFilter() const { return *this->filter; }
    const ScalingBasis &getScalingBasis() const { return this->basis; }
    const BoundingBox<D> &getWorldBox() const { return this->world; }

    MultiResolutionAnalysis<1> getKernelMRA() const;
    MultiResolutionAnalysis<2> getOperatorMRA() const;

    bool operator==(const MultiResolutionAnalysis<D> &mra) const;
    bool operator!=(const MultiResolutionAnalysis<D> &mra) const;

    void print() const;

protected:
    const int maxDepth;
    const ScalingBasis basis;
    const BoundingBox<D> world;
    MWFilter *filter;

    void setupFilter();
};

} // namespace mrcpp
