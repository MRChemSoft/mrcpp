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

#include "catch2/catch_all.hpp"

#include "MRCPP/constants.h"
#include "core/FilterCache.h"

using namespace Eigen;
using namespace mrcpp;

namespace mw_filter {

TEST_CASE("Interpolating filters", "[mw_filter]") {
    int maxOrder = 40;
    getInterpolatingFilterCache(ifilters);

    SECTION("Unit matrix") {
        double off_diag_row = 0.0;
        double off_diag_col = 0.0;
        for (int k = 1; k < maxOrder; k++) {
            for (int i = 0; i < k + 1; i++) {
                for (int j = i; j < k + 1; j++) {
                    const MatrixXd &F = ifilters.get(k).getFilter();
                    double sc = F.col(i).dot(F.col(j));
                    double sr = F.row(i).dot(F.row(j));
                    if (i == j) {
                        REQUIRE(sc == Catch::Approx(1.0));
                        REQUIRE(sr == Catch::Approx(1.0));
                    } else {
                        off_diag_col += std::abs(sc);
                        off_diag_row += std::abs(sr);
                    }
                }
            }
        }
        REQUIRE(off_diag_col == Catch::Approx(0.0).margin(1.0e-12));
        REQUIRE(off_diag_row == Catch::Approx(0.0).margin(1.0e-12));
    }
}

} // namespace mw_filter
