#include "catch.hpp"

#include "factory_functions.h"
#include "MWProjector.h"
#include "GridCleaner.h"
#include "GridGenerator.h"

namespace grid_cleaner {

template<int D> void testGridCleaner();

SCENARIO("Projected trees can be cleaned and reused", "[grid_cleaner], [tree_builder], [trees]") {
    GIVEN("A projected FunctionTree in 1D") {
        testGridCleaner<1>();
    }
    GIVEN("A projected FunctionTree in 2D") {
        testGridCleaner<2>();
    }
    GIVEN("A projected FunctionTree in 3D") {
        testGridCleaner<3>();
    }
}

template<int D> void testGridCleaner() {
    GaussFunc<D> *func = 0;
    initialize(&func);
    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    const double prec = 1.0e-4;
    GridGenerator<D> G;
    MWProjector<D> Q(prec);
    FunctionTree<D> tree(*mra);
    G(tree, *func);
    Q(tree, *func);

    const int refDepth = tree.getDepth();
    const int refNodes = tree.getNNodes();
    const double refInt = tree.integrate();
    const double refNorm = tree.getSquareNorm();

    WHEN("the grid is cleaned") {
        GridCleaner<D> C;
        C(tree);
        THEN("it represents an undefined function on the same grid") {
            REQUIRE( (tree.getDepth() == refDepth) );
            REQUIRE( (tree.getNNodes() == refNodes) );
            REQUIRE( (tree.integrate() == Approx(0.0)) );
            REQUIRE( (tree.getSquareNorm() == Approx(-1.0)) );
            AND_WHEN("the function is re-projected on the same grid") {
                Q(tree, *func);
                THEN("the representation is the same as before") {
                    REQUIRE( (tree.getDepth() == refDepth) );
                    REQUIRE( (tree.getNNodes() == refNodes) );
                    REQUIRE( (tree.integrate() == Approx(refInt)) );
                    REQUIRE( (tree.getSquareNorm() == Approx(refNorm)) );
                }
            }
        }
    }
    WHEN("the grid is cleaned adaptively") {
        const double new_prec = 1.0e-5;
        GridCleaner<D> C(new_prec);
        C(tree);
        THEN("it represents an undefined function on a larger grid") {
            REQUIRE( (tree.getDepth() >= refDepth) );
            REQUIRE( (tree.getNNodes() > refNodes) );
            REQUIRE( (tree.integrate() == Approx(0.0)) );
            REQUIRE( (tree.getSquareNorm() == Approx(-1.0)) );
            AND_WHEN("the function is re-projected on the new grid") {
                Q(tree, *func);
                THEN("it becomes a larger representation of the same function") {
                    REQUIRE( (tree.getDepth() >= refDepth) );
                    REQUIRE( (tree.getNNodes() > refNodes) );
                    REQUIRE( (tree.integrate() == Approx(refInt).epsilon(1.0e-8)) );
                    REQUIRE( (tree.getSquareNorm() == Approx(refNorm).epsilon(1.0e-8)) );
                }
            }
        }
    }
    finalize(&mra);
    finalize(&func);
}

} // namespace
