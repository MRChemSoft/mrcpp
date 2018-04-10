#include "catch.hpp"

#include "factory_functions.h"

#include "mwbuilders/grid.h"

using namespace mrcpp;

namespace grid_build {

template<int D> void testBuildGrid();

SCENARIO("Building empty grids", "[build_grid], [tree_builder], [trees]") {
    GIVEN("an analytic function in 1D") {
        testBuildGrid<1>();
    }
    GIVEN("an analytic function in 2D") {
        testBuildGrid<2>();
    }
    GIVEN("an analytic function in 3D") {
        testBuildGrid<3>();
    }
}

template<int D> void testBuildGrid() {
    GaussFunc<D> *f_func = 0;
    initialize(&f_func);

    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    WHEN("the GridGenerator is given no argument") {
        FunctionTree<D> f_tree(*mra);

        THEN("we get an empty tree of root nodes") {
            REQUIRE( (f_tree.getSquareNorm() == Approx(-1.0)) );
            REQUIRE( (f_tree.getDepth() == 1) );
            REQUIRE( (f_tree.getNNodes() == f_tree.getNEndNodes()) );
            REQUIRE( (f_tree.getNGenNodes() == 0) );

            AND_WHEN("the GridGenerator is given the analytic function") {
                build_grid(f_tree, *f_func, 2);

                THEN("the empty tree gets adapted") {
                    REQUIRE( (f_tree.getSquareNorm() == Approx(-1.0)) );
                    REQUIRE( (f_tree.getDepth() == 3) );
                    REQUIRE( (f_tree.getNNodes() > f_tree.getNEndNodes()) );
                    REQUIRE( (f_tree.getNGenNodes() == 0) );
                }
            }
        }
    }

    WHEN("the GridGenerator is given the analytic function") {
        FunctionTree<D> f_tree(*mra);
        build_grid(f_tree, *f_func, 2);

        THEN("we get an empty adapted tree structure") {
            REQUIRE( (f_tree.getSquareNorm() == Approx(-1.0)) );
            REQUIRE( (f_tree.getDepth() == 3) );
            REQUIRE( (f_tree.getNNodes() > f_tree.getNEndNodes()) );
            REQUIRE( (f_tree.getNGenNodes() == 0) );

            AND_WHEN("the empty tree is passed to the GridGenerator") {
                FunctionTree<D> g_tree(*mra);
                copy_grid(g_tree, f_tree);

                THEN("we get an identical empty grid") {
                    REQUIRE( (g_tree.getSquareNorm() == Approx(-1.0)) );
                    REQUIRE( (g_tree.getDepth() == f_tree.getDepth()) );
                    REQUIRE( (g_tree.getNNodes() == f_tree.getNNodes()) );
                    REQUIRE( (g_tree.getNEndNodes() == f_tree.getNEndNodes()) );
                    REQUIRE( (g_tree.getNGenNodes() == 0) );
                }
            }
        }
    }
    finalize(&mra);
    finalize(&f_func);
}

} // namespace
