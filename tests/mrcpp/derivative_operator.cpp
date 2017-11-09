#include "catch.hpp"

#include "factory_functions.h"


#include "MathUtils.h"
#include "ABGVOperator.h"
#include "MWAdder.h"
#include "MWDerivative.h"
#include "MWOperator.h"
#include "MWProjector.h"

using namespace std;

namespace derivative_operator {
template<int D>
MultiResolutionAnalysis<D>* initializeMRA() {
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
    return new MultiResolutionAnalysis<D>(world, basis);
}

template<int D> void testDifferentiation(double a,double b){

    MultiResolutionAnalysis<D> *mra = initializeMRA<D>();

    double prec = 1.0e-3;
    int max_scale = mra->getMaxScale();
    double proj_prec = prec/10.0;

    ABGVOperator<D> Differentiate(*mra, a, b);
    MWDerivative<D> apply(max_scale);

    MWAdder<D> add(-1.0, max_scale);
    MWProjector<D> project(proj_prec, max_scale);

    auto f = [] (const double *r) -> double {
        const double r_0[3] = {pi, pi, pi};
        double R = MathUtils::calcDistance(D, r, r_0);
        return exp(-R*R);
    };

    auto df = [] (const double *r) -> double {
        const double r_0[3] = {pi, pi, pi};
        double R = MathUtils::calcDistance(D, r, r_0);
        return -2*exp(-R*R)*(r[0]-pi);
    };

    FunctionTree<D> f_tree(*mra);
    project(f_tree, f);

    FunctionTree<D> df_tree(*mra);
    project(df_tree, df);

    FunctionTree<D> dg_tree(*mra);
    apply(dg_tree, Differentiate, f_tree, 0);

    FunctionTree<D> err_tree(*mra);
    add(err_tree, 1.0, df_tree, -1.0, dg_tree);

    double df_norm = sqrt(df_tree.getSquareNorm());
    double abs_err = sqrt(err_tree.getSquareNorm());
    double rel_err = abs_err/df_norm;

    REQUIRE(rel_err <= prec);


    delete mra;

}

TEST_CASE("differentiantion_entral_difference", "[derivative],[central_difference]") {
// 0.5,0.5 specifies central difference
    SECTION("1D_derivative_test"){
        testDifferentiation<1>(0.5,0.5);
    }
    SECTION("2D_derivative_test"){
        testDifferentiation<2>(0.5,0.5);
    }
    SECTION("3D derivative test"){
        testDifferentiation<3>(0.5,0.5);
    }
}
TEST_CASE("differentiantion_center_difference", "[derivative], [center_difference]") {
// 0,0 specifies center difference
    SECTION("1D_derivative_test"){
        testDifferentiation<1>(0,0);
    }
    SECTION("2D_derivative_test"){
        testDifferentiation<2>(0,0);
    }
    SECTION("3D_derivative_test"){
        testDifferentiation<3>(0,0);
    }
}
} // namespace
