#include "catch.hpp"

#include "factory_functions.h"
#include "PoissonOperator.h"
#include "MWOperator.h"
#include "MWConvolution.h"
#include "OperatorAdaptor.h"
#include "MWProjector.h"
#include "BandWidth.h"
#include "CrossCorrelationCalculator.h"
#include "GaussFunc.h"

using namespace std;

namespace poisson_operator {

TEST_CASE("Initialize Poisson operator", "[init_poisson], [poisson_operator], [mw_operator]") {
    const double r_min = 1.0e-2;
    const double r_max = 1.0e+1;
    const double exp_prec  = 1.0e-4;
    const double proj_prec = 1.0e-3;
    const double ccc_prec  = 1.0e-3;
    const double band_prec  = 1.0e-3;

    const int n = -3;
    const int k = 5;

    SECTION("Initialize Poisson's kernel") {
        PoissonKernel poisson(exp_prec, r_min, r_max);
        REQUIRE( (poisson.size() == 26) );

        double x = r_min;
        while (x < r_max) {
            REQUIRE( (poisson.evalf(&x) == Approx(1.0/x).epsilon(exp_prec)) );
            x *= 1.5;
        }
        SECTION("Project Poisson's kernel") {
            int l = -1;
            int nbox = 2;
            NodeIndex<1> idx(n, &l);
            BoundingBox<1> box(idx, &nbox);

            InterpolatingBasis basis(2*k+1);
            MultiResolutionAnalysis<1> kern_mra(box, basis);

            MWProjector<1> Q(proj_prec);
            GridGenerator<1> G;

            FunctionTreeVector<1> kern_vec;
            for (int i = 0; i < poisson.size(); i++) {
                Gaussian<1> &kern_gauss = *poisson[i];
                FunctionTree<1> *kern_tree = new FunctionTree<1>(kern_mra);
                G(*kern_tree, kern_gauss);
                Q(*kern_tree, kern_gauss);
                kern_vec.push_back(kern_tree);
            }

            SECTION("Build operator tree by cross correlation") {
                NodeIndex<2> idx(n);
                BoundingBox<2> box(idx);

                InterpolatingBasis basis(k);
                MultiResolutionAnalysis<2> oper_mra(box, basis);

                TreeBuilder<2> builder;
                OperatorAdaptor adaptor(ccc_prec, oper_mra.getMaxScale());

                MWOperator O(oper_mra);
                for (int i = 0; i < kern_vec.size(); i++) {
                    FunctionTree<1> &kern_tree = *kern_vec[i];
                    CrossCorrelationCalculator calculator(kern_tree);

                    OperatorTree *oper_tree = new OperatorTree(oper_mra, ccc_prec);
                    builder.build(*oper_tree, calculator, adaptor, -1);
                    oper_tree->setupOperNodeCache();
                    O.push_back(oper_tree);

                    oper_tree->calcBandWidth(1.0);
                    BandWidth bw_1 = oper_tree->getBandWidth();
                    oper_tree->clearBandWidth();

                    oper_tree->calcBandWidth(0.001);
                    BandWidth bw_2 = oper_tree->getBandWidth();
                    oper_tree->clearBandWidth();

                    oper_tree->calcBandWidth(-1.0);
                    BandWidth bw_3 = oper_tree->getBandWidth();
                    oper_tree->clearBandWidth();

                    for (int i = 0; i < oper_tree->getDepth(); i++) {
                        REQUIRE( (bw_1.getMaxWidth(i) <= bw_2.getMaxWidth(i)) );
                        REQUIRE( (bw_2.getMaxWidth(i) <= bw_3.getMaxWidth(i)) );
                    }
                }
                O.calcBandWidths(band_prec);
                REQUIRE( (O.getMaxBandWidth(3) == 3) );
                REQUIRE( (O.getMaxBandWidth(7) == 5) );
                REQUIRE( (O.getMaxBandWidth(13) == 9) );
                REQUIRE( (O.getMaxBandWidth(15) == -1) );

                O.clear(true);
            }
            for (int i = 0; i < kern_vec.size(); i++) {
                delete kern_vec[i];
            }
            kern_vec.clear();
        }
    }
}

TEST_CASE("Apply Poisson's operator", "[apply_poisson], [poisson_operator], [mw_operator]") {
    double proj_prec = 1.0e-4;
    double apply_prec = 1.0e-3;
    double build_prec = 1.0e-4;

    MultiResolutionAnalysis<3> *mra = 0;
    GaussFunc<3> *fFunc = 0;

    initialize(&fFunc);
    initialize(&mra);

    MWProjector<3> Q(proj_prec);
    PoissonOperator P(*mra, build_prec);
    MWConvolution<3> apply(apply_prec);

    FunctionTree<3> fTree(*mra);
    FunctionTree<3> gTree(*mra);

    Q(fTree, *fFunc);
    apply(gTree, P, fTree);

    double E_num = gTree.dot(fTree);
    double E_ana = fFunc->calcCoulombEnergy(*fFunc);

    REQUIRE( (E_num == Approx(E_ana).epsilon(apply_prec)) );

    finalize(&fFunc);
    finalize(&mra);
}

} // namespace
