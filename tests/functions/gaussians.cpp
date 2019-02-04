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

#include "MRCPP/functions/GaussPoly.h"
#include "MRCPP/treebuilders/grid.h"
#include "MRCPP/treebuilders/project.h"

using namespace mrcpp;

namespace gaussians {

template <int D> auto alpha_gen(const std::array<double, D> &sigma) {
    auto alpha = std::array<double, D>{};
    for (auto i = 0; i < D; i++) { alpha[i] = 0.5 / (sigma[i] * sigma[i]); }
    return alpha;
}

template <int D> auto coef_gen(const std::array<double, D> &sigma) {
    auto coef = 1.0;
    for (auto &s : sigma) coef *= 1.0 / (s * std::sqrt(2.0 * pi));
    return coef;
}

SCENARIO("Gaussians", "[gaussians]") {
    const auto D = 3;
    const auto prec = 1.0e-3;
    // Making normalized gaussian centered at the origin
    const auto sigma = std::array<double, D>{1.0, 2.0, 3.0};
    const auto pos = Coord<D>{};
    const auto power = std::array<int, D>{};
    const auto alpha = alpha_gen<D>(sigma);
    const auto coef = coef_gen<D>(sigma);
    auto gauss = GaussFunc<D>(alpha, coef, pos, power);

    // Making ref value
    const auto r_ref = std::array<double, D>{.1, .2, .3};
    const auto ref_val = coef * std::exp(-alpha[0] * r_ref[0] * r_ref[0] - alpha[1] * r_ref[1] * r_ref[1] -
                                         alpha[2] * r_ref[2] * r_ref[2]);

    // Making MRA
    const auto min_scale = -4;
    const auto corner = std::array<int, D>{-1, -1, -1};
    const auto boxes = std::array<int, D>{2, 2, 2};
    auto world = BoundingBox<D>(min_scale, corner, boxes);
    const auto basis = InterpolatingBasis(7);
    auto MRA = MultiResolutionAnalysis<D>(world, basis, 25);

    // Making a function tree and projects the gaussian onto it
    FunctionTree<D> f_tree(MRA);
    build_grid<D>(f_tree, gauss);
    project(prec, f_tree, gauss);

    WHEN("Gaussians are projected") {
        THEN("The gaussian can be evaluated") { REQUIRE(gauss.evalf(r_ref) == Approx(ref_val)); }

        THEN("the integral is normalized") { REQUIRE(f_tree.integrate() == Approx(1.0)); }

        THEN("multiply") {
            auto prod_gauss = gauss.mult(gauss);
            REQUIRE(prod_gauss.evalf(r_ref) == ref_val * ref_val);
        }
    }
}

} // namespace gaussians
