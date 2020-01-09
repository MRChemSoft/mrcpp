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

#include "catch.hpp"

#include "functions/LegendrePoly.h"

using namespace mrcpp;

namespace legendre_poly {

TEST_CASE("Legendre polynomials", "[legendre_poly], [polynomials]") {
    int nLeg = 10;
    std::vector<LegendrePoly *> L;
    for (int k = 0; k < nLeg; k++) {
        auto *L_k = new LegendrePoly(k, 2.0, 1.0);
        L.push_back(L_k);
    }

    SECTION("LegendrePoly constructor") {
        for (int k = 0; k < nLeg; k++) {
            LegendrePoly &L_k = *L[k];
            REQUIRE(L_k.getScaledLowerBound() == Approx(0.0));
            REQUIRE(L_k.getScaledUpperBound() == Approx(1.0));
            REQUIRE(L_k.getOrder() == k);
            // Legendre polynomials are normalized so that L_k(1.0) = 1.0
            Coord<1> one{1.0};
            REQUIRE(L_k.evalf(one) == Approx(1.0));
        }
    }

    SECTION("The Legendre polynomials form an orthogonal set") {
        for (int i = 0; i < nLeg; i++) {
            LegendrePoly &L_i = *L[i];
            double S_ii = L_i.innerProduct(L_i);
            REQUIRE(std::abs(S_ii) > MachineZero);
            for (int j = 0; j < i; j++) {
                LegendrePoly &L_j = *L[j];
                double S_ij = L_i.innerProduct(L_j);
                REQUIRE(S_ij == Approx(0.0).margin(1.0e-12));
            }
        }
    }
    for (int k = 0; k < nLeg; k++) { delete L[k]; }
}

} // namespace legendre_poly
