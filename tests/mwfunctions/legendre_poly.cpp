#include "catch.hpp"

#include "mwfunctions/LegendrePoly.h"

using namespace mrcpp;

namespace legendre_poly {

TEST_CASE("Legendre polynomials", "[legendre_poly], [polynomials]") {
    int nLeg = 10;
    std::vector<LegendrePoly *> L;
    for (int k = 0; k < nLeg; k++) {
        LegendrePoly *L_k = new LegendrePoly(k, 2.0, 1.0);
        L.push_back(L_k);
    }

    SECTION("LegendrePoly constructor") {
        for (int k = 0; k < nLeg; k++) {
            LegendrePoly &L_k = *L[k];
            REQUIRE( L_k.getScaledLowerBound() == Approx(0.0) );
            REQUIRE( L_k.getScaledUpperBound() == Approx(1.0) );
            REQUIRE( L_k.getOrder() == k );
            // Legendre polynomials are normalized so that L_k(1.0) = 1.0
            REQUIRE( L_k.evalf(1.0) == Approx(1.0) );
        }
    }

    SECTION("The Legendre polynomials form an orthogonal set") {
        for (int i = 0; i < nLeg; i++) {
            LegendrePoly &L_i = *L[i];
            double S_ii = L_i.innerProduct(L_i);
            REQUIRE( std::abs(S_ii) > MachineZero );
            for (int j = 0; j < i; j++) {
                LegendrePoly &L_j = *L[j];
                double S_ij = L_i.innerProduct(L_j);
                REQUIRE( S_ij == Approx(0.0).margin(1.0e-12) );
            }
        }
    }
    for (int k = 0; k < nLeg; k++) {
        delete L[k];
    }
}

} // namespace
