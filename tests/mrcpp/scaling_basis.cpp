#include "catch.hpp"

#include "InterpolatingBasis.h"
#include "LegendreBasis.h"

namespace scaling_basis {


template<int K> void testConstructor(const ScalingBasis &basis);
template<int K> void testOrthonormality(const ScalingBasis &basis);

TEST_CASE("InterpolatingBasis", "[interpolating_basis], [scaling_basis]") {
    SECTION("InterpolatingBasis order 1") {
        InterpolatingBasis basis(1);
        SECTION("Constructor") {
            for (int k = 0; k < 1; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE( (P_k.getScaledLowerBound() == Approx(0.0)) );
                REQUIRE( (P_k.getScaledUpperBound() == Approx(1.0)) );
                REQUIRE( (P_k.getOrder() == 1) );
            }
        }
        SECTION("Orthonormality") {
            testOrthonormality<1>(basis);
        }
    }
    SECTION("InterpolatingBasis order 6") {
        InterpolatingBasis basis(6);
        SECTION("Constructor") {
            for (int k = 0; k < 6; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE( (P_k.getScaledLowerBound() == Approx(0.0)) );
                REQUIRE( (P_k.getScaledUpperBound() == Approx(1.0)) );
                REQUIRE( (P_k.getOrder() == 6) );
            }
        }
        SECTION("Orthonormality") {
            testOrthonormality<6>(basis);
        }
    }
    SECTION("InterpolatingBasis order 10") {
        InterpolatingBasis basis(10);
        SECTION("Constructor") {
            for (int k = 0; k < 10; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE( (P_k.getScaledLowerBound() == Approx(0.0)) );
                REQUIRE( (P_k.getScaledUpperBound() == Approx(1.0)) );
                REQUIRE( (P_k.getOrder() == 10) );
            }
        }
        SECTION("Orthonormality") {
            testOrthonormality<10>(basis);
        }
    }
}

TEST_CASE("LegendreBasis", "[legendre_basis], [scaling_basis]") {
    SECTION("LegendreBasis order 1") {
        LegendreBasis basis(1);
        SECTION("Constructor") {
            for (int k = 0; k < 1; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE( (P_k.getScaledLowerBound() == Approx(0.0)) );
                REQUIRE( (P_k.getScaledUpperBound() == Approx(1.0)) );
                REQUIRE( (P_k.getOrder() == k) );
            }
        }
        SECTION("Orthonormality") {
            testOrthonormality<1>(basis);
        }
    }
    SECTION("LegendreBasis order 6") {
        LegendreBasis basis(6);
        SECTION("Test constructor") {
            for (int k = 0; k < 6; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE( (P_k.getScaledLowerBound() == Approx(0.0)) );
                REQUIRE( (P_k.getScaledUpperBound() == Approx(1.0)) );
                REQUIRE( (P_k.getOrder() == k) );
            }
        }
        SECTION("Orthonormality") {
            testOrthonormality<6>(basis);
        }
    }
    SECTION("LegendreBasis order 10") {
        LegendreBasis basis(10);
        SECTION("Test constructor") {
            for (int k = 0; k < 10; k++) {
                const Polynomial &P_k = basis.getFunc(k);
                REQUIRE( (P_k.getScaledLowerBound() == Approx(0.0)) );
                REQUIRE( (P_k.getScaledUpperBound() == Approx(1.0)) );
                REQUIRE( (P_k.getOrder() == k) );
            }
        }
        SECTION("Orthonormality") {
            testOrthonormality<10>(basis);
        }
    }
}

template<int K>
void testOrthonormality(const ScalingBasis &basis) {
    for (int i = 0; i < K; i++) {
        const Polynomial &P_i = basis.getFunc(i);
        double S_ii = P_i.innerProduct(P_i);
        REQUIRE( (S_ii == Approx(1.0)) );
        for (int j = 0; j < i; j++) {
            const Polynomial &P_j = basis.getFunc(j);
            double S_ij = P_i.innerProduct(P_j);
            REQUIRE( (fabs(S_ij) == Approx(0.0)) );
        }
    }
}
} // namespace
