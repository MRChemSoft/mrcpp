#include "catch.hpp"

#include "factory_functions.h"

#include "TreeBuilder.h"
#include "HelmholtzOperator.h"
#include "HelmholtzKernel.h"
#include "MWOperator.h"
#include "OperatorAdaptor.h"
#include "BandWidth.h"
#include "CrossCorrelationCalculator.h"
#include "HydrogenicFunction.h"
#include "MathUtils.h"
#include "apply.h"
#include "project.h"
#include "multiply.h"
#include "add.h"
#include "grid.h"

using namespace std;
using namespace mrcpp;

namespace helmholtz_operator {

TEST_CASE("Helmholtz' kernel", "[init_helmholtz], [helmholtz_operator], [mw_operator]") {
    const double mu = 1.0;
    const double r_min = 1.0e-3;
    const double r_max = 1.0e+1;
    const double exp_prec  = 1.0e-4;
    const double proj_prec = 1.0e-3;
    const double ccc_prec  = 1.0e-3;
    const double band_prec  = 1.0e-3;

    const int n = -3;
    const int k = 5;

    SECTION("Initialize Helmholtz' kernel") {
        HelmholtzKernel helmholtz(mu, exp_prec, r_min, r_max);
        REQUIRE( (helmholtz.size() == 24) );

        int foo = 0;
        double x = r_min;
        while (x < r_max) {
            REQUIRE( (helmholtz.evalf(&x) == Approx(exp(-mu*x)/x).epsilon(exp_prec)) );
            x *= 1.5;
        }
        SECTION("Project Helmholtz' kernel") {
            int l = -1;
            int nbox = 2;
            NodeIndex<1> idx(n, &l);
            BoundingBox<1> box(idx, &nbox);

            InterpolatingBasis basis(2*k+1);
            MultiResolutionAnalysis<1> kern_mra(box, basis);

            FunctionTreeVector<1> K;
            for (int i = 0; i < helmholtz.size(); i++) {
                Gaussian<1> &kern_gauss = *helmholtz[i];
                FunctionTree<1> *kern_tree = new FunctionTree<1>(kern_mra);
                build_grid(*kern_tree, kern_gauss);
                project(proj_prec, *kern_tree, kern_gauss);
                K.push_back(kern_tree);
            }

            SECTION("Build operator tree by cross correlation") {
                NodeIndex<2> idx(n);
                BoundingBox<2> box(idx);

                InterpolatingBasis basis(k);
                MultiResolutionAnalysis<2> oper_mra(box, basis);

                TreeBuilder<2> builder;
                OperatorAdaptor adaptor(ccc_prec, oper_mra.getMaxScale());

                MWOperator O(oper_mra);
                for (int i = 0; i < K.size(); i++) {
                    FunctionTree<1> &kern_tree = *K[i];
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
                REQUIRE( (O.getMaxBandWidth(20) == -1) );

                O.clear(true);
            }
            K.clear(true);
        }
    }
}

TEST_CASE("Apply Helmholtz' operator", "[apply_helmholtz], [helmholtz_operator], [mw_operator]") {
    double proj_prec = 3.0e-3;
    double apply_prec = 3.0e-2;
    double build_prec = 3.0e-3;

    // Computational domain [-32.0, 32.0]
    int scale = -5;
    int corner[3] = {-1, -1, -1};
    int nbox[3] = {2, 2, 2};
    NodeIndex<3> idx(scale, corner);
    BoundingBox<3> box(idx, nbox);

    int order = 5;
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(box, basis);

    int n = 2;                  // Principal quantum number
    int l = 1;                  // Angular quantum number
    int m_l = 2;                // Magnetic quantum number
    double Z = 1.0;             // Nuclear charge
    double E = -Z/(2.0*n*n);    // Total energy

    double mu = sqrt(-2*E);
    HelmholtzOperator H(MRA, mu, build_prec);

    double R[3] = {0.0, 0.0, 0.0};
    HydrogenicFunction hFunc(n, l, m_l, Z, R);
    FunctionTree<3> psi_n(MRA);
    project(proj_prec, psi_n, hFunc);

    auto f = [Z, R] (const double *r) -> double {
        double x = MathUtils::calcDistance(3, r, R);
        return -Z/x;
    };
    FunctionTree<3> V(MRA);
    project(proj_prec, V, f);

    FunctionTree<3> Vpsi(MRA);
    copy_grid(Vpsi, psi_n);
    multiply(-1.0, Vpsi, 1.0, V, psi_n);

    FunctionTree<3> psi_np1(MRA);
    copy_grid(psi_np1, psi_n);
    apply(apply_prec, psi_np1, H, Vpsi);
    psi_np1.rescale(-1.0/(2.0*pi));

    double norm = sqrt(psi_np1.getSquareNorm());
    REQUIRE( (norm == Approx(1.0).epsilon(apply_prec)) );

    FunctionTree<3> d_psi(MRA);
    copy_grid(d_psi, psi_np1);
    add(-1.0, d_psi, 1.0, psi_np1, -1.0, psi_n);

    double error = sqrt(d_psi.getSquareNorm());
    REQUIRE( (error < apply_prec) );
}

} // namespace
