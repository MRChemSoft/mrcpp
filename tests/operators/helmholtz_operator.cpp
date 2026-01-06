/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "catch2/catch_all.hpp"

#include "factory_functions.h"

#include "operators/HelmholtzKernel.h"
#include "operators/HelmholtzOperator.h"
#include "operators/MWOperator.h"
#include "treebuilders//OperatorAdaptor.h"
#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/add.h"
#include "treebuilders/apply.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"
#include "trees/BandWidth.h"
#include "utils/math_utils.h"

using namespace mrcpp;

namespace helmholtz_operator {

TEST_CASE("Helmholtz' kernel", "[init_helmholtz], [helmholtz_operator], [mw_operator]") {
    const double mu = 0.01;
    const double r_min = 1.0e-3;
    const double r_max = 1.0e+0;
    const double exp_prec = 1.0e-4;
    const double proj_prec = 1.0e-3;
    const double ccc_prec = 1.0e-3;

    const int n = -3;
    const int k = 5;

    NodeIndex<3> func_idx(n);
    BoundingBox<3> func_box(func_idx);

    InterpolatingBasis func_basis(k);
    MultiResolutionAnalysis<3> func_mra(func_box, func_basis);

    SECTION("Initialize Helmholtz' kernel") {
        HelmholtzKernel helmholtz(mu, exp_prec, r_min, r_max);
        REQUIRE(helmholtz.size() == 33);

        Coord<1> x{r_min};
        while (x[0] < r_max) {
            REQUIRE(helmholtz.evalf(x) == Catch::Approx(std::exp(-mu * x[0]) / x[0]).epsilon(2.0 * exp_prec));
            x[0] *= 1.5;
        }
        SECTION("Project Helmholtz' kernel") {
            NodeIndex<1> kern_idx(n, {-1});
            BoundingBox<1> kern_box(kern_idx, {2});

            InterpolatingBasis kern_basis(2 * k + 1);
            MultiResolutionAnalysis<1> kern_mra(kern_box, kern_basis);

            FunctionTreeVector<1> K;
            for (int i = 0; i < helmholtz.size(); i++) {
                Gaussian<1> &kern_gauss = *helmholtz[i];
                auto *kern_tree = new FunctionTree<1>(kern_mra);
                build_grid(*kern_tree, kern_gauss);
                project(proj_prec, *kern_tree, kern_gauss);
                K.push_back(std::make_tuple(1.0, kern_tree));
            }

            SECTION("Build operator tree by cross correlation") {
                NodeIndex<2> oper_idx(n);
                BoundingBox<2> oper_box(oper_idx);

                InterpolatingBasis oper_basis(k);
                MultiResolutionAnalysis<2> oper_mra(oper_box, oper_basis);

                TreeBuilder<2> builder;
                OperatorAdaptor adaptor(ccc_prec, oper_mra.getMaxScale());

                for (int i = 0; i < K.size(); i++) {
                    FunctionTree<1> &kern_tree = get_func(K, i);
                    CrossCorrelationCalculator calculator(kern_tree);

                    auto oper_tree = std::make_unique<OperatorTree>(oper_mra, ccc_prec);
                    builder.build(*oper_tree, calculator, adaptor, -1);
                    oper_tree->setupOperNodeCache();

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
                        REQUIRE(bw_1.getMaxWidth(i) <= bw_2.getMaxWidth(i));
                        REQUIRE(bw_2.getMaxWidth(i) <= bw_3.getMaxWidth(i));
                    }
                }
            }
            clear(K, true);
        }
    }
}

TEST_CASE("Apply Helmholtz' operator", "[apply_helmholtz], [helmholtz_operator], [mw_operator]") {
    double proj_prec = 3.0e-3;
    double apply_prec = 3.0e-2;
    double build_prec = 3.0e-3;

    // Computational domain [-32.0, 32.0]
    int scale = -4;
    std::array<int, 3> corner;
    std::array<int, 3> nbox;
    std::array<double, 3> scaling_factor;
    corner.fill(-1);
    nbox.fill(2);
    scaling_factor.fill(2.0);
    BoundingBox<3> box(scale, corner, nbox, scaling_factor);
    int order = 5;

    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(box, basis);

    int n = 1;                     // Principal quantum number
    double Z = 1.0;                // Nuclear charge
    double E = -Z / (2.0 * n * n); // Total energy

    double mu = std::sqrt(-2 * E);

    double r_min = MRA.calcMinDistance(build_prec);
    double r_max = MRA.calcMaxDistance();

    HelmholtzKernel helmholtz(mu, build_prec, r_min, r_max);
    ConvolutionOperator<3> H(MRA, helmholtz, build_prec);

    // Defining analytic 1s function
    auto hFunc = [Z](const Coord<3> &r) -> double {
        const double c_0 = 2.0 * std::pow(Z, 3.0 / 2.0);
        double rho = 2.0 * Z * std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        double R_0 = c_0 * std::exp(-rho / 2.0);
        double Y_00 = 1.0 / std::sqrt(4.0 * mrcpp::pi);
        return R_0 * Y_00;
    };
    FunctionTree<3> psi_n(MRA);
    project<3, double>(proj_prec, psi_n, hFunc);

    auto f = [Z](const Coord<3> &r) -> double {
        double x = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return -Z / x;
    };
    FunctionTree<3> V(MRA);
    project<3, double>(proj_prec, V, f);

    FunctionTree<3> Vpsi(MRA);
    copy_grid(Vpsi, psi_n);
    multiply(-1.0, Vpsi, 1.0, V, psi_n);

    FunctionTree<3> psi_np1(MRA);
    copy_grid(psi_np1, psi_n);
    apply(apply_prec, psi_np1, H, Vpsi);
    psi_np1.rescale(-1.0 / (2.0 * pi));

    double norm = std::sqrt(psi_np1.getSquareNorm());
    REQUIRE(norm == Catch::Approx(1.0).epsilon(apply_prec));

    FunctionTree<3> d_psi(MRA);
    copy_grid(d_psi, psi_np1);
    add(-1.0, d_psi, 1.0, psi_np1, -1.0, psi_n);

    double error = std::sqrt(d_psi.getSquareNorm());
    REQUIRE(error == Catch::Approx(0.0).margin(apply_prec));
}

TEST_CASE("Apply Periodic Helmholtz' operator", "[apply_periodic_helmholtz], [helmholtz_operator], [mw_operator]") {
    double proj_prec = 3.0e-3;
    double apply_prec = 3.0e-2;
    double build_prec = 3.0e-3;

    // 2.0*pi periodic in all dirs // UPDATE ME
    auto scaling_factor = std::array<double, 3>{pi, pi, pi};
    auto corner = std::array<int, 3>{-1, -1, -1};
    auto boxes = std::array<int, 3>{2, 2, 2};
    auto world = mrcpp::BoundingBox<3>(0, corner, boxes, scaling_factor, true);
    int order = 5;

    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(world, basis, 25);

    auto mu = 4.3;
    auto oper_root = 0;
    auto oper_reach = 9;
    HelmholtzOperator H(MRA, mu, build_prec, oper_root, oper_reach);

    // Source, Poisson applied to this should yield cos(x)cos(y)cos(z)
    auto source = [mu](const mrcpp::Coord<3> &r) { return 3.0 * cos(r[0]) * cos(r[1]) * cos(r[2]) / (4.0 * pi) + mu * mu * cos(r[0]) * cos(r[1]) * cos(r[2]) / (4.0 * pi); };

    FunctionTree<3> source_tree(MRA);
    project<3, double>(proj_prec, source_tree, source);

    FunctionTree<3> sol_tree(MRA);
    FunctionTree<3> in_tree(MRA);
    FunctionTree<3> out_tree(MRA);
    FunctionTree<3> in_out_tree(MRA);

    apply(apply_prec, sol_tree, H, source_tree);
    apply_near_field(apply_prec, in_tree, H, source_tree);
    apply_far_field(apply_prec, out_tree, H, source_tree);

    add(apply_prec, in_out_tree, 1.0, in_tree, 1.0, out_tree);

    REQUIRE(sol_tree.evalf({0.0, 0.0, 0.0}) == Catch::Approx(1.0).epsilon(apply_prec));
    REQUIRE(sol_tree.evalf({pi, 0.0, 0.0}) == Catch::Approx(-1.0).epsilon(apply_prec));
    REQUIRE(in_out_tree.evalf({0.0, 0.0, 0.0}) == Catch::Approx(1.0).epsilon(apply_prec));
    REQUIRE(in_out_tree.evalf({pi, 0.0, 0.0}) == Catch::Approx(-1.0).epsilon(apply_prec));
}

TEST_CASE("Apply negative scale Helmholtz' operator", "[apply_periodic_helmholtz], [helmholtz_operator], [mw_operator]. [negative_scale]") {
    double proj_prec = 3.0e-3;
    double apply_prec = 3.0e-2;
    double build_prec = 3.0e-3;

    // 2.0*pi periodic in all dirs // UPDATE ME
    auto scaling_factor = std::array<double, 3>{pi, pi, pi};
    auto corner = std::array<int, 3>{-1, -1, -1};
    auto boxes = std::array<int, 3>{2, 2, 2};
    auto world = mrcpp::BoundingBox<3>(0, corner, boxes, scaling_factor, true);
    int order = 5;

    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(world, basis, 25);

    auto mu = 4.3;
    auto oper_root = -4;
    auto oper_cutoff = 1;
    HelmholtzOperator H(MRA, mu, build_prec, oper_root, oper_cutoff);

    // Source, Poisson applied to this should yield cos(x)cos(y)cos(z)
    auto source = [mu](const mrcpp::Coord<3> &r) { return 3.0 * cos(r[0]) * cos(r[1]) * cos(r[2]) / (4.0 * pi) + mu * mu * cos(r[0]) * cos(r[1]) * cos(r[2]) / (4.0 * pi); };

    FunctionTree<3> source_tree(MRA);
    project<3, double>(proj_prec, source_tree, source);

    FunctionTree<3> sol_tree(MRA);

    apply(apply_prec, sol_tree, H, source_tree);

    REQUIRE(sol_tree.evalf({0.0, 0.0, 0.0}) == Catch::Approx(1.0).epsilon(apply_prec));
    REQUIRE(sol_tree.evalf({pi, 0.0, 0.0}) == Catch::Approx(-1.0).epsilon(apply_prec));
}

} // namespace helmholtz_operator
