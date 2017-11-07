#include "catch.hpp"

#include "parallel.h"
#include "Timer.h"
#include "MathUtils.h"
#include "ABGVOperator.h"
#include "MWAdder.h"
#include "MWDerivative.h"

#include "factory_functions.h"
#include "MWOperator.h"
#include "MWConvolution.h"
#include "OperatorAdaptor.h"
#include "MWProjector.h"
#include "BandWidth.h"
#include "CrossCorrelationCalculator.h"
#include "GaussFunc.h"

using namespace std;

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
    return new MultiResolutionAnalysis<D>(world, basis);
}

TEST_CASE("Differentiantion", "[derivative]") {



    SECTION("1D derivative test") {

        MultiResolutionAnalysis<1> *MRA_1 = initializeMRA<1>();

        double prec = 1.0e-3;
        int max_scale = MRA_1->getMaxScale();
        double proj_prec = prec/10.0;

        ABGVOperator<1> D(*MRA_1, 0.0, 0.0);
        MWDerivative<1> apply(max_scale);

        MWAdder<1> add(-1.0, max_scale);
        MWProjector<1> project(proj_prec, max_scale);

        auto f = [] (const double *r) -> double {
            const double alpha = 3.0;
            const double r_0[3] = {pi, pi, pi};
            double R = MathUtils::calcDistance(1, r, r_0);
            return exp(-alpha*R);
        };
        auto df = [] (const double *r) -> double {
            const double alpha = 3.0;
            const double r_0[3] = {pi, pi, pi};
            double R = MathUtils::calcDistance(1, r, r_0);
            double sign = 1.0;
            if (r[0] > r_0[0]) sign = -1.0;
            return sign*alpha*exp(-alpha*R);
        };

        FunctionTree<1> f_tree(*MRA_1);
        project(f_tree, f);

        FunctionTree<1> df_tree(*MRA_1);
        project(df_tree, df);

        FunctionTree<1> dg_tree(*MRA_1);
        apply(dg_tree, D, f_tree, 0); // Does not refine grid further

        FunctionTree<1> err_tree(*MRA_1);
        add(err_tree, 1.0, df_tree, -1.0, dg_tree);

        double df_norm = sqrt(df_tree.getSquareNorm());
        double abs_err = sqrt(err_tree.getSquareNorm());
        double rel_err = abs_err/df_norm;

        delete MRA_1;


        REQUIRE(rel_err<=prec);
    }
}

} // namespace
