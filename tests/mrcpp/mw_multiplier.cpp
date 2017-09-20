#include "catch.hpp"

#include "factory_functions.h"
#include "GridGenerator.h"
#include "MWProjector.h"
#include "MWMultiplier.h"
#include "WaveletAdaptor.h"
#include "GaussPoly.h"

namespace mw_multiplier {

template<int D> void testMultiplication();

SCENARIO("MWMultiplier", "[mw_multiplier], [tree_builder]") {
    GIVEN("Two MW functions in 1D") {
        testMultiplication<1>();
    }
    GIVEN("Two MW functions in 2D") {
        testMultiplication<2>();
    }
    GIVEN("Two MW functions in 3D") {
        testMultiplication<3>();
    }
}

template<int D> void testMultiplication() {
    double alpha = 1.0;
    double beta_a = 110.0;
    double beta_b = 50.0;
    double pos_a[3] = {-0.25, 0.35, 1.05};
    double pos_b[3] = {-0.20, 0.50, 1.05};

    GaussFunc<D> a_func(beta_a, alpha, pos_a);
    GaussFunc<D> b_func(beta_b, alpha, pos_b);
    GaussPoly<D> ref_func = a_func*b_func;

    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    // Setting up adaptor and TreeBuilders
    double prec = 1.0e-4;
    GridGenerator<D> G;
    MWProjector<D> Q(prec);
    MWMultiplier<D> mult;

    // Initialize trees
    FunctionTree<D> a_tree(*mra);
    FunctionTree<D> b_tree(*mra);
    FunctionTree<D> ref_tree(*mra);

    // Build empty grids
    G(a_tree, a_func);
    G(b_tree, b_func);
    G(ref_tree, ref_func);

    // Project functions
    Q(a_tree, a_func);
    Q(b_tree, b_func);
    Q(ref_tree, ref_func);

    const double ref_int = ref_tree.integrate();
    const double ref_norm = ref_tree.getSquareNorm();

    FunctionTreeVector<D> prod_vec;
    WHEN("the functions are multiplied") {
        FunctionTree<D> c_tree(*mra);
        prod_vec.push_back(&a_tree);
        prod_vec.push_back(&b_tree);
        mult(c_tree, prod_vec);
        prod_vec.clear();

        THEN("the MW product equals the analytic product") {
            double c_int = c_tree.integrate();
            double c_dot = c_tree.dot(ref_tree);
            double c_norm = c_tree.getSquareNorm();
            REQUIRE( (c_int == Approx(ref_int)) );
            REQUIRE( (c_dot == Approx(ref_norm)) );
            REQUIRE( (c_norm == Approx(ref_norm)) );
        }
    }
    finalize(&mra);
}

} // namespace
