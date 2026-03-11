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

#include "functions/GaussFunc.h"
#include "operators/MWOperator.h"
#include "operators/PoissonKernel.h"
#include "operators/PoissonOperator.h"
#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/apply.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"
#include "trees/BandWidth.h"

using namespace mrcpp;

namespace poisson_operator {

TEST_CASE("Initialize Poisson operator", "[init_poisson], [poisson_operator], [mw_operator]") {
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

    SECTION("Initialize Poisson's kernel") {
        PoissonKernel poisson(exp_prec, r_min, r_max);
        REQUIRE(poisson.size() == 26);

        Coord<1> x{r_min};
        while (x[0] < r_max) {
            REQUIRE(poisson.evalf(x) == Catch::Approx(1.0 / x[0]).epsilon(2.0 * exp_prec));
            x[0] *= 1.5;
        }
        SECTION("Project Poisson's kernel") {
            NodeIndex<1> kern_idx(n, {-1});
            BoundingBox<1> kern_box(kern_idx, {2});

            InterpolatingBasis kern_basis(2 * k + 1);
            MultiResolutionAnalysis<1> kern_mra(kern_box, kern_basis);

            FunctionTreeVector<1> kern_vec;
            for (int i = 0; i < poisson.size(); i++) {
                Gaussian<1> &kern_gauss = *poisson[i];
                auto *kern_tree = new FunctionTree<1>(kern_mra);
                build_grid(*kern_tree, kern_gauss);
                project(proj_prec, *kern_tree, kern_gauss);
                kern_vec.push_back(std::make_tuple(1.0, kern_tree));
            }

            SECTION("Build operator tree by cross correlation") {
                NodeIndex<2> idx(n);
                BoundingBox<2> box(idx);

                InterpolatingBasis basis(k);
                MultiResolutionAnalysis<2> oper_mra(box, basis);

                TreeBuilder<2> builder;
                OperatorAdaptor adaptor(ccc_prec, oper_mra.getMaxScale());

                for (int i = 0; i < kern_vec.size(); i++) {
                    FunctionTree<1> &kern_tree = get_func(kern_vec, i);
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
            clear(kern_vec, true);
        }
    }
}

TEST_CASE("Apply Poisson's operator", "[apply_poisson], [poisson_operator], [mw_operator]") {
    double proj_prec = 1.0e-4;
    double apply_prec = 1.0e-3;
    double build_prec = 1.0e-4;

    MultiResolutionAnalysis<3> *mra = nullptr;
    GaussFunc<3> *fFunc = nullptr;

    initialize(&fFunc);
    initialize(&mra);

    PoissonOperator P(*mra, build_prec);
    FunctionTree<3> fTree(*mra);
    FunctionTree<3> gTree(*mra);

    project(proj_prec, fTree, *fFunc);
    apply(apply_prec, gTree, P, fTree);

    double E_num = dot(gTree, fTree);
    double E_ana = fFunc->calcCoulombEnergy(*fFunc);

    REQUIRE(E_num == Catch::Approx(E_ana).epsilon(apply_prec));

    finalize(&fFunc);
    finalize(&mra);
}

TEST_CASE("Apply Periodic Poisson' operator", "[apply_periodic_Poisson], [poisson_operator], [mw_operator]") {
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

    auto oper_scale = 0;
    auto oper_reach = 9;
    auto r_min = MRA.calcMinDistance(build_prec);
    auto r_max = MRA.calcMaxDistance();

    // Regular non-periodic operators should have rel_root=0 and oper_reach=0
    auto rel_root = oper_scale - MRA.getRootScale();
    r_max *= std::pow(2.0, -rel_root);
    r_max *= (2.0 * oper_reach) + 1.0;

    PoissonKernel poisson(build_prec, r_min, r_max);
    ConvolutionOperator<3> P(MRA, poisson, build_prec, oper_scale, oper_reach);

    // Source, Poisson applied to this should yield cos(x)cos(y)cos(z)
    auto source = [](const mrcpp::Coord<3> &r) { return 3.0 * cos(r[0]) * cos(r[1]) * cos(r[2]) / (4.0 * pi); };

    FunctionTree<3> source_tree(MRA);
    project<3, double>(proj_prec, source_tree, source);

    FunctionTree<3> sol_tree(MRA);

    apply(apply_prec, sol_tree, P, source_tree);

    REQUIRE(sol_tree.evalf({0.0, 0.0, 0.0}) == Catch::Approx(1.0).epsilon(apply_prec));
    REQUIRE(sol_tree.evalf({pi, 0.0, 0.0}) == Catch::Approx(-1.0).epsilon(apply_prec));
}

} // namespace poisson_operator
