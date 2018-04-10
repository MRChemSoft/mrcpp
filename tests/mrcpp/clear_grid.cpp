#include "catch.hpp"

#include "factory_functions.h"

#include "mwbuilders/project.h"
#include "mwbuilders/grid.h"

using namespace mrcpp;

namespace grid_clear {

template<int D> void testClearGrid();

SCENARIO("Projected trees can be cleared and reused", "[clear_grid], [tree_builder], [trees]") {
    GIVEN("A projected FunctionTree in 1D") {
        testClearGrid<1>();
    }
    GIVEN("A projected FunctionTree in 2D") {
        testClearGrid<2>();
    }
    GIVEN("A projected FunctionTree in 3D") {
        testClearGrid<3>();
    }
}

template<int D> void testClearGrid() {
    GaussFunc<D> *func = 0;
    initialize(&func);
    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    const double prec = 1.0e-4;
    FunctionTree<D> tree(*mra);
    build_grid(tree, *func);
    project(prec, tree, *func);

    const int refDepth = tree.getDepth();
    const int refNodes = tree.getNNodes();
    const double refInt = tree.integrate();
    const double refNorm = tree.getSquareNorm();

    WHEN("the tree is cleared") {
        tree.clear();
        THEN("it represents an undefined function on the root grid") {
            REQUIRE( (tree.getDepth() == 1) );
            REQUIRE( (tree.getNNodes() == tree.getRootBox().size()) );
            REQUIRE( (tree.integrate() == Approx(0.0)) );
            REQUIRE( (tree.getSquareNorm() == Approx(-1.0)) );
            AND_WHEN("the function is re-projected") {
                build_grid(tree, *func);
                project(prec, tree, *func);
                THEN("the representation is the same as before it was cleared") {
                    REQUIRE( (tree.getDepth() == refDepth) );
                    REQUIRE( (tree.getNNodes() == refNodes) );
                    REQUIRE( (tree.integrate() == Approx(refInt)) );
                    REQUIRE( (tree.getSquareNorm() == Approx(refNorm)) );
                }
            }
        }
    }
    WHEN("the grid is cleared") {
        clear_grid(-1.0, tree);
        THEN("it represents an undefined function on the same grid") {
            REQUIRE( (tree.getDepth() == refDepth) );
            REQUIRE( (tree.getNNodes() == refNodes) );
            REQUIRE( (tree.integrate() == Approx(0.0)) );
            REQUIRE( (tree.getSquareNorm() == Approx(-1.0)) );
            AND_WHEN("the function is re-projected on the same grid") {
                project(-1.0, tree, *func);
                THEN("the representation is the same as before") {
                    REQUIRE( (tree.getDepth() == refDepth) );
                    REQUIRE( (tree.getNNodes() == refNodes) );
                    REQUIRE( (tree.integrate() == Approx(refInt)) );
                    REQUIRE( (tree.getSquareNorm() == Approx(refNorm)) );
                }
            }
        }
    }
    WHEN("the grid is cleared adaptively") {
        const double new_prec = 1.0e-5;
        clear_grid(new_prec, tree);
        THEN("it represents an undefined function on a larger grid") {
            REQUIRE( (tree.getDepth() >= refDepth) );
            REQUIRE( (tree.getNNodes() > refNodes) );
            REQUIRE( (tree.integrate() == Approx(0.0)) );
            REQUIRE( (tree.getSquareNorm() == Approx(-1.0)) );
            AND_WHEN("the function is re-projected on the new grid") {
                project(-1.0, tree, *func);
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
