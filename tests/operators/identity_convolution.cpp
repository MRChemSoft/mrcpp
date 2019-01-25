/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "catch.hpp"

#include "factory_functions.h"

#include "operators/IdentityConvolution.h"
#include "operators/IdentityKernel.h"
#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/apply.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"
#include "trees/BandWidth.h"
#include "trees/OperatorTree.h"

using namespace mrcpp;

namespace identity_convolution {

template <int D> void applyIdentity();

TEST_CASE("Initialize identity convolution operator", "[init_identity], [identity_convolution], [mw_operator]") {
    double exp_prec = 1.0e-6;
    double proj_prec = 1.0e-6;
    double ccc_prec = 1.0e-4;

    const int n = -6;
    const int k = 5;

    SECTION("Initialize identity kernel") {
        IdentityKernel id_kern(exp_prec);
        REQUIRE(id_kern.size() == 1);

        SECTION("Project identity kernel") {
            std::array<int, 1> l{-1};
            std::array<int, 1> nbox{2};
            NodeIndex<1> idx(n, l.data());
            BoundingBox<1> box(idx, nbox);

            InterpolatingBasis basis(2 * k + 1);
            MultiResolutionAnalysis<1> kern_mra(box, basis);

            FunctionTree<1> kern_tree(kern_mra);
            build_grid(kern_tree, id_kern);
            project(proj_prec, kern_tree, id_kern);
            REQUIRE(kern_tree.integrate() == Approx(1.0).epsilon(proj_prec));

            SECTION("Build operator tree by cross correlation") {
                NodeIndex<2> idx(n);
                BoundingBox<2> box(idx);

                InterpolatingBasis basis(k);
                MultiResolutionAnalysis<2> oper_mra(box, basis);

                TreeBuilder<2> builder;
                OperatorAdaptor adaptor(ccc_prec, oper_mra.getMaxScale());
                CrossCorrelationCalculator calculator(kern_tree);

                OperatorTree oper_tree(oper_mra, ccc_prec);
                builder.build(oper_tree, calculator, adaptor, -1);
                oper_tree.setupOperNodeCache();

                oper_tree.calcBandWidth(1.0);
                BandWidth bw_1 = oper_tree.getBandWidth();
                oper_tree.clearBandWidth();

                oper_tree.calcBandWidth(0.001);
                BandWidth bw_2 = oper_tree.getBandWidth();
                oper_tree.clearBandWidth();

                oper_tree.calcBandWidth(-1.0);
                BandWidth bw_3 = oper_tree.getBandWidth();
                oper_tree.clearBandWidth();

                for (int i = 0; i < oper_tree.getDepth(); i++) {
                    REQUIRE(bw_1.getMaxWidth(i) <= bw_2.getMaxWidth(i));
                    REQUIRE(bw_2.getMaxWidth(i) <= bw_3.getMaxWidth(i));
                }
            }
        }
    }
}

TEST_CASE("Apply identity convolution operator", "[apply_identity], [identity_convolution], [mw_operator]") {
    SECTION("1D") { applyIdentity<1>(); }
    SECTION("2D") { applyIdentity<2>(); }
    SECTION("3D") { applyIdentity<3>(); }
}

template <int D> void applyIdentity() {
    double proj_prec = 1.0e-3;
    double apply_prec = 1.0e-3;
    double build_prec = 1.0e-4;

    MultiResolutionAnalysis<D> *mra = 0;
    GaussFunc<D> *fFunc = 0;
    initialize(&fFunc);
    initialize(&mra);

    IdentityConvolution<D> I(*mra, build_prec);
    FunctionTree<D> fTree(*mra);
    FunctionTree<D> gTree(*mra);

    project(proj_prec, fTree, *fFunc);
    apply(apply_prec, gTree, I, fTree);

    REQUIRE(gTree.getDepth() <= fTree.getDepth());
    REQUIRE(gTree.getNNodes() <= fTree.getNNodes());
    REQUIRE(gTree.integrate() == Approx(fTree.integrate()).epsilon(apply_prec));

    finalize(&fFunc);
    finalize(&mra);
}

} // namespace identity_convolution
