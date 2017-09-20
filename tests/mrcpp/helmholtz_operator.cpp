#include "catch.hpp"

#include "factory_functions.h"
#include "HelmholtzOperator.h"
#include "MWOperator.h"
#include "MWConvolution.h"
#include "OperatorAdaptor.h"
#include "MWProjector.h"
#include "MWMultiplier.h"
#include "MWAdder.h"
#include "BandWidth.h"
#include "CrossCorrelationCalculator.h"
#include "HydrogenicFunction.h"
#include "MathUtils.h"

using namespace std;

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

            MWProjector<1> Q(proj_prec);
            GridGenerator<1> G;

            FunctionTreeVector<1> K;
            for (int i = 0; i < helmholtz.size(); i++) {
                Gaussian<1> &kern_gauss = *helmholtz[i];
                FunctionTree<1> *kern_tree = new FunctionTree<1>(kern_mra);
                G(*kern_tree, kern_gauss);
                Q(*kern_tree, kern_gauss);
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
    double proj_prec = 1.0e-3;
    double apply_prec = 1.0e-2;
    double build_prec = 1.0e-3;

    // Computational domain [-32.0, 32.0]
    int scale = -5;
    int corner[3] = {-1, -1, -1};
    int nbox[3] = {2, 2, 2};
    NodeIndex<3> idx(scale, corner);
    BoundingBox<3> box(idx, nbox);

    int order = 5;
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(box, basis);

    FunctionTreeVector<3> tree_vec;
    MWAdder<3> add;
    MWMultiplier<3> mult;
    MWProjector<3> Q(proj_prec);
    GridGenerator<3> G;
    MWConvolution<3> apply(apply_prec);

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
    Q(psi_n, hFunc);

    auto f = [Z, R] (const double *r) -> double {
        double x = MathUtils::calcDistance(3, r, R);
        return -Z/x;
    };
    FunctionTree<3> V(MRA);
    Q(V, f);

    FunctionTree<3> Vpsi(MRA);
    G(Vpsi, psi_n);
    mult(Vpsi, 1.0, V, psi_n);

    FunctionTree<3> psi_np1(MRA);
    G(psi_np1, psi_n);
    apply(psi_np1, H, Vpsi);
    psi_np1 *= -1.0/(2.0*pi);

    double norm = sqrt(psi_np1.getSquareNorm());
    REQUIRE( (norm == Approx(1.0).epsilon(apply_prec)) );

    FunctionTree<3> d_psi(MRA);
    G(d_psi, psi_np1);
    add(d_psi, 1.0, psi_np1, -1.0, psi_n);

    double error = sqrt(d_psi.getSquareNorm());
    REQUIRE( (error < apply_prec) );
}

} // namespace
