#include "catch.hpp"

#include "factory_functions.h"
#include "GridGenerator.h"
#include "MWProjector.h"
#include "MWAdder.h"
#include "WaveletAdaptor.h"
#include "GaussExp.h"

namespace mw_adder {

template<int D> void testAddition();

SCENARIO("Adding MW trees", "[mw_adder], [tree_builder]") {
    GIVEN("Two MW functions in 1D") {
        testAddition<1>();
    }
    GIVEN("Two MW functions in 2D") {
        testAddition<2>();
    }
    GIVEN("Two MW functions in 3D") {
        testAddition<3>();
    }
}

template<int D> void testAddition() {
    const double a_coef = 1.0;
    const double b_coef = 2.0;

    double alpha = 1.0;
    double beta_a = 110.0;
    double beta_b = 50.0;
    double pos_a[3] = {-0.25, 0.35, 1.05};
    double pos_b[3] = {-0.20, 0.50, 1.05};

    GaussFunc<D> a_func(beta_a, alpha, pos_a);
    GaussFunc<D> b_func(beta_b, alpha, pos_b);
    GaussExp<D> ref_func;
    ref_func.append(a_func);
    ref_func.append(b_func);
    ref_func.setCoef(0, a_coef);
    ref_func.setCoef(1, b_coef);

    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    // Setting up adaptor and TreeBuilders
    double prec = 1.0e-4;
    GridGenerator<D> G;
    MWProjector<D> Q(prec);
    MWAdder<D> add;

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

    // Reference integrals
    const double a_int = a_tree.integrate();
    const double b_int = b_tree.integrate();

    const double ref_int = ref_tree.integrate();
    const double ref_norm = ref_tree.getSquareNorm();

    FunctionTreeVector<D> sum_vec;
    WHEN("the functions are added") {
        FunctionTree<D> c_tree(*mra);
        sum_vec.push_back(a_coef, &a_tree);
        sum_vec.push_back(b_coef, &b_tree);
        add(c_tree, sum_vec);
        sum_vec.clear();

        THEN("their integrals add up") {
            double c_int = c_tree.integrate();
            double int_sum = a_coef*a_int + b_coef*b_int;
            REQUIRE( (c_int == Approx(int_sum)) );
        }

        AND_THEN("the MW sum equals the analytic sum") {
            double c_int = c_tree.integrate();
            double c_dot = c_tree.dot(ref_tree);
            double c_norm = c_tree.getSquareNorm();
            REQUIRE( (c_int == Approx(ref_int)) );
            REQUIRE( (c_dot == Approx(ref_norm)) );
            REQUIRE( (c_norm == Approx(ref_norm)) );
        }

        AND_WHEN("the first function is subtracted") {
            FunctionTree<D> d_tree(*mra);
            sum_vec.push_back(&c_tree);
            sum_vec.push_back(-1.0, &a_tree);
            add(d_tree, sum_vec);
            sum_vec.clear();

            THEN("the integral is the same as the second function") {
                double d_int = d_tree.integrate();
                double ref_int = b_coef*b_int;
                REQUIRE( (d_int == Approx(ref_int)) );
            }

            AND_WHEN("the second function is subtracted") {
                FunctionTree<D> e_tree(*mra);
                sum_vec.push_back(&d_tree);
                sum_vec.push_back(-b_coef, &b_tree);
                add(e_tree, sum_vec);
                sum_vec.clear();

                THEN("the integral is zero") {
                    double e_int = e_tree.integrate();
                    double ref_int = 0.0;
                    REQUIRE( (e_int == Approx(ref_int)) );
                }
            }
        }
    }
    finalize(&mra);
}

} // namespace
