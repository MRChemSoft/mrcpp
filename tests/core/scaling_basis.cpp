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

#include "catch.hpp"

#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"

using namespace mrcpp;

namespace scaling_basis {

template <int K> void testConstructor(const ScalingBasis &basis);
template <int K> void testOrthonormality(const ScalingBasis &basis);

TEST_CASE("InterpolatingBasis", "[interpolating_basis], [scaling_basis]") {
    SECTION("InterpolatingBasis order 1") {
        InterpolatingBasis basis(1);
        SECTION("Constructor") {
            for (int k = 0; k < 1; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE(P_k.getScaledLowerBound() == Approx(0.0));
                REQUIRE(P_k.getScaledUpperBound() == Approx(1.0));
                REQUIRE(P_k.getOrder() == 1);
            }
        }
        SECTION("Orthonormality") { testOrthonormality<1>(basis); }
    }
    SECTION("InterpolatingBasis order 6") {
        InterpolatingBasis basis(6);
        SECTION("Constructor") {
            for (int k = 0; k < 6; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE(P_k.getScaledLowerBound() == Approx(0.0));
                REQUIRE(P_k.getScaledUpperBound() == Approx(1.0));
                REQUIRE(P_k.getOrder() == 6);
            }
        }
        SECTION("Orthonormality") { testOrthonormality<6>(basis); }
    }
    SECTION("InterpolatingBasis order 10") {
        InterpolatingBasis basis(10);
        SECTION("Constructor") {
            for (int k = 0; k < 10; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE(P_k.getScaledLowerBound() == Approx(0.0));
                REQUIRE(P_k.getScaledUpperBound() == Approx(1.0));
                REQUIRE(P_k.getOrder() == 10);
            }
        }
        SECTION("Orthonormality") { testOrthonormality<10>(basis); }
    }
}

TEST_CASE("LegendreBasis", "[legendre_basis], [scaling_basis]") {
    SECTION("LegendreBasis order 1") {
        LegendreBasis basis(1);
        SECTION("Constructor") {
            for (int k = 0; k < 1; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE(P_k.getScaledLowerBound() == Approx(0.0));
                REQUIRE(P_k.getScaledUpperBound() == Approx(1.0));
                REQUIRE(P_k.getOrder() == k);
            }
        }
        SECTION("Orthonormality") { testOrthonormality<1>(basis); }
    }
    SECTION("LegendreBasis order 6") {
        LegendreBasis basis(6);
        SECTION("Test constructor") {
            for (int k = 0; k < 6; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE(P_k.getScaledLowerBound() == Approx(0.0));
                REQUIRE(P_k.getScaledUpperBound() == Approx(1.0));
                REQUIRE(P_k.getOrder() == k);
            }
        }
        SECTION("Orthonormality") { testOrthonormality<6>(basis); }
    }
    SECTION("LegendreBasis order 10") {
        LegendreBasis basis(10);
        SECTION("Test constructor") {
            for (int k = 0; k < 10; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE(P_k.getScaledLowerBound() == Approx(0.0));
                REQUIRE(P_k.getScaledUpperBound() == Approx(1.0));
                REQUIRE(P_k.getOrder() == k);
            }
        }
        SECTION("Orthonormality") { testOrthonormality<10>(basis); }
    }
}

template <int K> void testOrthonormality(const ScalingBasis &basis) {
    for (int i = 0; i < K; i++) {
        const Polynomial &P_i = basis.getFunc(i);
        double S_ii = P_i.innerProduct(P_i);
        REQUIRE(S_ii == Approx(1.0));
        for (int j = 0; j < i; j++) {
            const Polynomial &P_j = basis.getFunc(j);
            double S_ij = P_i.innerProduct(P_j);
            REQUIRE(S_ij == Approx(0.0).margin(1.0e-10));
        }
    }
}
} // namespace scaling_basis
