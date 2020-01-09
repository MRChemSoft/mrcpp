/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "treebuilders/grid.h"
#include "treebuilders/project.h"

using namespace mrcpp;

namespace scaling {

template <int D> void testScaling();

template <typename T, int D> constexpr auto generate_array(const T &inp) {
    auto arr = std::array<T, D>{};
    arr.fill(inp);
    return arr;
}

template <int D> auto alpha_gen(double sigma) {
    return std::pow(1.0 / std::sqrt(2.0 * pi * std::pow(sigma, 2.0)), D);
}

auto beta_gen(double sigma) {
    return 1.0 / (2.0 * std::pow(sigma, 2.0));
}

TEST_CASE("Scaling factor integral", "[scaling_factor]") {
    SECTION("1D") { testScaling<1>(); }
    SECTION("2D") { testScaling<2>(); }
    SECTION("3D") { testScaling<3>(); }
}

template <int D> void testScaling() {
    const auto prec = 1.0e-3;
    const auto min_scale = 0;

    const auto corner = std::array<int, D>{};
    const auto boxes = generate_array<int, D>(1);
    const auto sf = generate_array<double, D>(2.0 * pi);

    const auto basis = InterpolatingBasis(5);
    const auto world = BoundingBox<D>(min_scale, corner, boxes, sf);
    const auto MRA = MultiResolutionAnalysis<D>(world, basis, 25);

    const auto sigma = 0.01;
    const auto alpha = alpha_gen<D>(sigma);
    const auto beta = beta_gen(sigma);
    const auto power = std::array<int, D>{};
    const auto pos = generate_array<double, D>(3.0);
    auto gauss = GaussFunc<D>(beta, alpha, pos, power);

    FunctionTree<D> f_tree(MRA);
    build_grid(f_tree, gauss);
    project(prec, f_tree, gauss);

    REQUIRE(Approx(1.0) == f_tree.integrate());
}

} // namespace scaling
