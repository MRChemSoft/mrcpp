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

/** Builds an OperatorTree with known band width (e.g. derivative and identity).
 * Assumes translational invariant and symmetric (in x - y) operator and keeps
 * only lower row of nodes (lx = 0). */

namespace mrcpp {

class BandWidthAdaptor final : public TreeAdaptor<2> {
public:
    BandWidthAdaptor(int bw, int ms)
            : TreeAdaptor<2>(ms)
            , bandWidth(bw) {}

private:
    const int bandWidth;

    bool splitNode(const MWNode<2> &node) const override {
        const auto &idx = node.getNodeIndex();
        int dl = std::abs(idx[0] - idx[1]);
        // Within band width on NEXT scale
        return ((idx[0] == 0) and (2 * dl <= this->bandWidth));
    }
};

} // namespace mrcpp
