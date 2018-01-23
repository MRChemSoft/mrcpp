#include "catch.hpp"

#include "factory_functions.h"

#include "FunctionTreeVector.h"
#include "WaveletAdaptor.h"
#include "GaussPoly.h"
#include "project.h"
#include "grid.h"
#include "multiply.h"

using namespace mrcpp;

namespace multiplication {

template<int D> void testMultiplication();

SCENARIO("Multiplying MW trees", "[multiplication], [tree_builder]") {
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
    const double prec = 1.0e-4;

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

    // Initialize trees
    FunctionTree<D> a_tree(*mra);
    FunctionTree<D> b_tree(*mra);
    FunctionTree<D> ref_tree(*mra);

    // Build empty grids
    build_grid(a_tree, a_func);
    build_grid(b_tree, b_func);
    build_grid(ref_tree, ref_func);

    // Project functions
    project(prec, a_tree, a_func);
    project(prec, b_tree, b_func);
    project(prec, ref_tree, ref_func);

    const double ref_int = ref_tree.integrate();
    const double ref_norm = ref_tree.getSquareNorm();

    FunctionTreeVector<D> prod_vec;
    WHEN("the functions are multiplied") {
        FunctionTree<D> c_tree(*mra);
        prod_vec.push_back(&a_tree);
        prod_vec.push_back(&b_tree);
        multiply(prec, c_tree, prod_vec);
        prod_vec.clear();

        THEN("the MW product equals the analytic product") {
            double c_int = c_tree.integrate();
            double c_dot = dot(c_tree, ref_tree);
            double c_norm = c_tree.getSquareNorm();
            REQUIRE( (c_int == Approx(ref_int)) );
            REQUIRE( (c_dot == Approx(ref_norm)) );
            REQUIRE( (c_norm == Approx(ref_norm)) );
        }
    }
    finalize(&mra);
}

} // namespace
