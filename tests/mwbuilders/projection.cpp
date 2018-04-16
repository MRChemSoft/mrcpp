#include "catch.hpp"

#include "factory_functions.h"

#include "mwbuilders/project.h"
#include "mwbuilders/grid.h"

using namespace mrcpp;

namespace projection {

template<int D> void testProjectFunction();

SCENARIO("Projecting Gaussian function", "[projection], [tree_builder], [trees]") {
    GIVEN("a Gaussian of unit charge in 1D") {
        testProjectFunction<1>();
    }
    GIVEN("a Gaussian of unit charge in 2D") {
        testProjectFunction<2>();
    }
    GIVEN("a Gaussian of unit charge in 3D") {
        testProjectFunction<3>();
    }
}

template<int D> void testProjectFunction() {
    GaussFunc<D> *func = 0;
    initialize(&func);
    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    WHEN("the function is projected on the default grid") {
        FunctionTree<D> tree(*mra);
        project(-1.0, tree, *func);
        THEN("it integrates to approximately one") {
            REQUIRE( tree.integrate() == Approx(1.0).margin(1.0) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = tree.getSquareNorm();
            REQUIRE( dot(tree, tree) == Approx(norm) );
        }
    }
    WHEN("the function is projected on an adapted grid") {
        FunctionTree<D> tree(*mra);
        build_grid(tree, *func);
        project(-1.0, tree, *func);
        THEN("it integrates to approximately one") {
            REQUIRE( tree.integrate() == Approx(1.0).epsilon(1.0e-3) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = tree.getSquareNorm();
            REQUIRE( dot(tree, tree) == Approx(norm) );
        }
    }
    WHEN("the function is projected with guaranteed precision") {
        const double prec = 1.0e-4;
        FunctionTree<D> f_tree(*mra);
        project(prec, f_tree, *func);
        THEN("it integrates to approximately one") {
            REQUIRE( f_tree.integrate() == Approx(1.0).epsilon(1.0e-8) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = f_tree.getSquareNorm();
            REQUIRE( dot(f_tree, f_tree) == Approx(norm) );
        }
        AND_WHEN("the function is projected on an identical grid") {
            FunctionTree<D> g_tree(*mra);
            copy_grid(g_tree, f_tree);

            project(-1.0, g_tree, *func);
            THEN("it integrates to the same value") {
                const double charge = f_tree.integrate();
                REQUIRE( g_tree.integrate() == Approx(charge) );
            }
            THEN("the dot product with the original is equal to their squared norm") {
                const double norm = f_tree.getSquareNorm();
                REQUIRE( dot(g_tree, f_tree) == Approx(norm) );
            }
        }
    }
    finalize(&mra);
    finalize(&func);
}

} // namespace
