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

#include "catch.hpp"

#include "trees/BandWidth.h"

using namespace Eigen;
using namespace mrcpp;

namespace band_width {

TEST_CASE("BandWidth", "[band_width]") {
    const int depth = 10;
    BandWidth bw(depth);

    SECTION("Default BandWidth") {
        REQUIRE(bw.getDepth() == depth);
        for (int n = 0; n < 2 * depth; n++) {
            REQUIRE(bw.isEmpty(n));
            REQUIRE(bw.getMaxWidth(n) == -1);
        }
    }

    SECTION("Setting bandnd widths") {
        bw.setWidth(0, 0, 1);
        bw.setWidth(0, 1, 2);
        bw.setWidth(0, 2, 3);
        bw.setWidth(0, 3, 0);

        bw.setWidth(1, 0, 5);
        bw.setWidth(1, 1, 4);
        bw.setWidth(1, 2, 3);
        bw.setWidth(1, 3, 1);

        REQUIRE(bw.getDepth() == depth);
        REQUIRE_FALSE(bw.isEmpty(0));
        REQUIRE_FALSE(bw.isEmpty(1));
        REQUIRE(bw.getMaxWidth(0) == 3);
        REQUIRE(bw.getMaxWidth(1) == 5);
        for (int n = 2; n < 2 * depth; n++) {
            REQUIRE(bw.isEmpty(n));
            REQUIRE(bw.getMaxWidth(n) == -1);
        }
    }
}

} // namespace band_width
