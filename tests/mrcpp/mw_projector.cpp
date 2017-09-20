#include "catch.hpp"

#include "factory_functions.h"
#include "MWProjector.h"

namespace mw_projector {

template<int D> void testProjectFunction();

SCENARIO("Projecting Gaussian function", "[mw_projector], [tree_builder], [trees]") {
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
        MWProjector<D> Q;
        FunctionTree<D> tree(*mra);
        Q(tree, *func);
        THEN("it integrates to approximately one") {
            REQUIRE( (tree.integrate() == Approx(1.0).epsilon(1.0e+1)) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = tree.getSquareNorm();
            REQUIRE( (tree.dot(tree) == Approx(norm)) );
        }
    }
    WHEN("the function is projected on an adapted grid") {
        GridGenerator<D> G;
        MWProjector<D> Q;
        FunctionTree<D> tree(*mra);
        G(tree, *func);
        Q(tree, *func);
        THEN("it integrates to approximately one") {
            REQUIRE( (tree.integrate() == Approx(1.0).epsilon(1.0e-3)) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = tree.getSquareNorm();
            REQUIRE( (tree.dot(tree) == Approx(norm)) );
        }
    }
    WHEN("the function is projected with guaranteed precision") {
        const double prec = 1.0e-4;
        MWProjector<D> Q_adap(prec);
        FunctionTree<D> f_tree(*mra);
        Q_adap(f_tree, *func);
        THEN("it integrates to approximately one") {
            REQUIRE( (f_tree.integrate() == Approx(1.0).epsilon(1.0e-8)) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = f_tree.getSquareNorm();
            REQUIRE( (f_tree.dot(f_tree) == Approx(norm)) );
        }
        AND_WHEN("the function is projected on an identical grid") {
            GridGenerator<D> G;
            FunctionTree<D> g_tree(*mra);
            G(g_tree, f_tree);

            MWProjector<D> Q_copy;
            Q_copy(g_tree, *func);
            THEN("it integrates to the same value") {
                const double charge = f_tree.integrate();
                REQUIRE( (g_tree.integrate() == Approx(charge)) );
            }
            THEN("the dot product with the original is equal to their squared norm") {
                const double norm = f_tree.getSquareNorm();
                REQUIRE( (g_tree.dot(f_tree) == Approx(norm)) );
            }
        }
    }
    finalize(&mra);
    finalize(&func);
}

} // namespace
