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
 * <https://mrcpp.readthedocs.io/en/latest/>
 */

#include "catch2/catch_all.hpp"

#include "factory_functions.h"

#include "operators/TimeEvolutionOperator.h"
#include "trees/MWNode.h"
#include "trees/OperatorTree.h"
#include "utils/math_utils.h"

#include <algorithm>
#include <cmath>
#include <complex>

using namespace mrcpp;

namespace complex_schrodinger_evolution_operator {

/** @brief The native complex Schrodinger kernel must equal the real cos + i*sin split.
 *
 *  The classic construction builds the Schrodinger semigroup exp(i t d^2) as two
 *  *real* operator trees: imaginary = false gives the cos part, imaginary = true
 *  the sin part. The new TimeEvolutionOperator<1, ComplexDouble> instead stores a
 *  single *complex* kernel whose coefficients are cos + i*sin.
 *
 *  Because the uniform tree builder fixes the operator-tree structure from the
 *  scale alone (independently of the coefficient values), the three operator trees
 *  are structurally identical, so we can compare them leaf by leaf. Every step that
 *  produces the coefficients (the cross-correlation contraction, the per-node MW
 *  compression, and the bottom-up MW transform) is linear with real filters, so the
 *  complex coefficients must equal (real coef) + i*(imag coef) to round-off.
 */
TEST_CASE("Complex Schrodinger evolution operator equals the real/imag split",
          "[complex_schrodinger_evolution_operator], [mw_operator]") {
    const auto min_scale = 0;
    const auto max_depth = 25;
    const auto order = 4;
    const auto prec = 1.0e-7;

    const int finest_scale = 7;
    const double delta_t = 0.03;

    auto basis = LegendreBasis(order);
    auto world = BoundingBox<1>(min_scale);
    auto MRA = MultiResolutionAnalysis<1>(world, basis, max_depth);

    // Classic split: two real operator trees (cos and sin parts).
    TimeEvolutionOperator<1, double> ReExp(MRA, prec, delta_t, finest_scale, false);
    TimeEvolutionOperator<1, double> ImExp(MRA, prec, delta_t, finest_scale, true);

    // Native single complex kernel (the bool flag is ignored for the complex type).
    TimeEvolutionOperator<1, ComplexDouble> Exp(MRA, prec, delta_t, finest_scale, false);

    const OperatorTree<double> &re = ReExp.getComponent(0, 0);
    const OperatorTree<double> &im = ImExp.getComponent(0, 0);
    const OperatorTree<ComplexDouble> &cx = Exp.getComponent(0, 0);

    // (1) Square-norm identity: ||cos + i sin||^2 == ||cos||^2 + ||sin||^2.
    const double sq_re = re.getSquareNorm();
    const double sq_im = im.getSquareNorm();
    const double sq_cx = cx.getSquareNorm();
    REQUIRE(sq_cx == Catch::Approx(sq_re + sq_im).epsilon(1.0e-10));

    // (2) Leaf-by-leaf: complex coef == real coef + i * imag coef.
    REQUIRE(cx.getNEndNodes() == re.getNEndNodes());
    REQUIRE(cx.getNEndNodes() == im.getNEndNodes());

    double max_diff = 0.0;
    for (int n = 0; n < cx.getNEndNodes(); n++) {
        const MWNode<2, ComplexDouble> &cx_node = cx.getEndMWNode(n);
        const MWNode<2, double> &re_node = re.getEndMWNode(n);
        const MWNode<2, double> &im_node = im.getEndMWNode(n);

        // The uniform builder must produce identical structure in all three trees.
        REQUIRE(cx_node.getNodeIndex() == re_node.getNodeIndex());
        REQUIRE(cx_node.getNodeIndex() == im_node.getNodeIndex());

        const int n_coefs = cx_node.getNCoefs();
        REQUIRE(re_node.getNCoefs() == n_coefs);
        REQUIRE(im_node.getNCoefs() == n_coefs);

        const ComplexDouble *cx_c = cx_node.getCoefs();
        const double *re_c = re_node.getCoefs();
        const double *im_c = im_node.getCoefs();
        for (int k = 0; k < n_coefs; k++) {
            const ComplexDouble expected(re_c[k], im_c[k]);
            max_diff = std::max(max_diff, std::abs(cx_c[k] - expected));
        }
    }

    REQUIRE(max_diff == Catch::Approx(0.0).margin(1.0e-12));
}

} // namespace complex_schrodinger_evolution_operator
