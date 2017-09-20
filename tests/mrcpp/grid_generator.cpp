#include "catch.hpp"

#include "factory_functions.h"

namespace grid_generator {

template<int D> void testGridGenerator();

SCENARIO("The GridGenerator builds empty grids", "[grid_generator], [tree_builder], [trees]") {
    GIVEN("a GridGenerator and analytic function in 1D") {
        testGridGenerator<1>();
    }
    GIVEN("a GridGenerator and analytic function in 2D") {
        testGridGenerator<2>();
    }
    GIVEN("a GridGenerator and analytic function in 3D") {
        testGridGenerator<3>();
    }
}

template<int D> void testGridGenerator() {
    GaussFunc<D> *f_func = 0;
    initialize(&f_func);

    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);
    GridGenerator<D> G;

    WHEN("the GridGenerator is given no argument") {
        FunctionTree<D> f_tree(*mra);

        THEN("we get an empty tree of root nodes") {
            REQUIRE( (f_tree.getSquareNorm() == Approx(-1.0)) );
            REQUIRE( (f_tree.getDepth() == 1) );
            REQUIRE( (f_tree.getNNodes() == f_tree.getNEndNodes()) );
            REQUIRE( (f_tree.getNGenNodes() == 0) );

            AND_WHEN("the GridGenerator is given the analytic function") {
                G(f_tree, *f_func, 2);

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
        G(f_tree, *f_func, 2);

        THEN("we get an empty adapted tree structure") {
            REQUIRE( (f_tree.getSquareNorm() == Approx(-1.0)) );
            REQUIRE( (f_tree.getDepth() == 3) );
            REQUIRE( (f_tree.getNNodes() > f_tree.getNEndNodes()) );
            REQUIRE( (f_tree.getNGenNodes() == 0) );

            AND_WHEN("the empty tree is passed to the GridGenerator") {
                FunctionTree<D> g_tree(*mra);
                G(g_tree, f_tree);

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
