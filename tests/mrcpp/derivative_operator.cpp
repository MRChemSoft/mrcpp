#include "catch.hpp"

#include "factory_functions.h"

#include "mwutils/MathUtils.h"
#include "mwoperators/MWOperator.h"
#include "mwoperators/ABGVOperator.h"
#include "mwoperators/PHOperator.h"
#include "mwbuilders/add.h"
#include "mwbuilders/apply.h"
#include "mwbuilders/project.h"

using namespace std;
using namespace mrcpp;

namespace derivative_operator {

template<int D> MultiResolutionAnalysis<D>* initializeMRA() {
    // Constructing world box
    int min_scale = -4;
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    NodeIndex<D> idx(min_scale, corner);
    BoundingBox<D> world(idx, boxes);

    // Constructing scaling basis
    int order = 5;
    InterpolatingBasis basis(order);

    // Initializing MRA
    int max_depth = 20;
    return new MultiResolutionAnalysis<D>(world, basis, max_depth);
}

template<int D> void testDifferentiationABGV(double a, double b) {
    MultiResolutionAnalysis<D> *mra = initializeMRA<D>();

    double prec = 1.0e-3;
    ABGVOperator<D> diff(*mra, a, b);

    const double r_0[3] = {pi, pi, pi};
    auto f = [r_0] (const double *r) -> double {
        double R = MathUtils::calcDistance(D, r, r_0);
        return exp(-R*R);
    };

    auto df = [r_0] (const double *r) -> double {
        double R = MathUtils::calcDistance(D, r, r_0);
        return -2.0*exp(-R*R)*(r[0]-r_0[0]);
    };

    FunctionTree<D> f_tree(*mra);
    project(prec/10, f_tree, f);

    FunctionTree<D> df_tree(*mra);
    project(prec/10, df_tree, df);

    FunctionTree<D> dg_tree(*mra);
    apply(dg_tree, diff, f_tree, 0);

    FunctionTree<D> err_tree(*mra);
    add(-1.0, err_tree, 1.0, df_tree, -1.0, dg_tree);

    double df_norm = sqrt(df_tree.getSquareNorm());
    double abs_err = sqrt(err_tree.getSquareNorm());
    double rel_err = abs_err/df_norm;

    REQUIRE(rel_err <= prec);

    delete mra;
}

template<int D> void testDifferentiationPH(int order) {
    MultiResolutionAnalysis<D> *mra = initializeMRA<D>();

    double prec = 1.0e-3;
    PHOperator<D> diff(*mra, order);

    const double r_0[3] = {pi, pi, pi};
    auto f = [r_0] (const double *r) -> double {
        double R = MathUtils::calcDistance(D, r, r_0);
        return exp(-R*R);
    };
    auto df = [r_0, order] (const double *r) -> double {
        double R = MathUtils::calcDistance(D, r, r_0);
        return -(2-order)*2*exp(-R*R)*(r[0]-pi)
                + (order-1)*(-2*exp(-R*R)+4*exp(-R*R)*(r[0]-r_0[0])*(r[0]-r_0[0]));
                // 2-order = 1 and order-1 = 0 in the first order case
                // 2-order = 0 and order-1 = 1 in the second order case
    };

    FunctionTree<D> f_tree(*mra);
    project(prec/10, f_tree, f);

    FunctionTree<D> df_tree(*mra);
    project(prec/10, df_tree, df);

    FunctionTree<D> dg_tree(*mra);
    apply(dg_tree, diff, f_tree, 0);

    FunctionTree<D> err_tree(*mra);
    add(-1.0, err_tree, 1.0, df_tree, -1.0, dg_tree);

    double df_norm = sqrt(df_tree.getSquareNorm());
    double abs_err = sqrt(err_tree.getSquareNorm());
    double rel_err = abs_err/df_norm;

    REQUIRE(rel_err <= prec);

    delete mra;
}

TEST_CASE("ABGV differentiantion central difference", "[derivative], [central_difference]") {
    // 0.5,0.5 specifies central difference
    SECTION("1D derivative test"){
        testDifferentiationABGV<1>(0.5, 0.5);
    }
    SECTION("2D derivative test"){
        testDifferentiationABGV<2>(0.5, 0.5);
    }
    SECTION("3D derivative test"){
        testDifferentiationABGV<3>(0.5, 0.5);
    }
}

TEST_CASE("ABGV differentiantion center difference", "[derivative], [center_difference]") {
    // 0,0 specifies center difference
    SECTION("1D derivative test") {
        testDifferentiationABGV<1>(0, 0);
    }
    SECTION("2D derivative test") {
        testDifferentiationABGV<2>(0, 0);
    }
    SECTION("3D derivative test") {
        testDifferentiationABGV<3>(0, 0);
    }
}

TEST_CASE("PH differentiantion first order", "[derivative], [PH_first_order]") {
    SECTION("1D derivative test") {
        testDifferentiationPH<1>(1);
    }
    SECTION("2D derivative test") {
        testDifferentiationPH<2>(1);
    }
    SECTION("3D derivative test") {
        testDifferentiationPH<3>(1);
    }
}

TEST_CASE("PH differentiantion second order", "[derivative], [PH_second_order]") {
    SECTION("1D second order derivative test") {
        testDifferentiationPH<1>(2);
    }
    SECTION("2D second order derivative test") {
        testDifferentiationPH<2>(2);
    }
    SECTION("3D second order derivative test") {
        testDifferentiationPH<3>(2);
    }
}

} // namespace
