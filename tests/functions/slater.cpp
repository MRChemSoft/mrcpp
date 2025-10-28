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

#include "functions/Slater.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"
#include "utils/math_utils.h"

using namespace mrcpp;

namespace gaussians {

template <int D> void testSlaterFunction();

SCENARIO("Slater", "[slater]") {
    GIVEN("Slater Function in 1D") { testSlaterFunction<1>(); }
    GIVEN("Slater Function in 2D") { testSlaterFunction<2>(); }
    GIVEN("Slater Function in 3D") { testSlaterFunction<3>(); }
}

template <int D> void testSlaterFunction() {
    const auto prec = 1.0e-4;
    // Making normalized gaussian centered at the origin
    const auto alpha = 2.;
    const auto coeff = 3.;   
    double pos_data[3] = {0.1, 0.2, 0.3};
    double ref_data[3] = {-0.2, 0.5, 1.0};
    auto pos = mrcpp::details::convert_to_std_array<double, D>(pos_data);
    auto ref = mrcpp::details::convert_to_std_array<double, D>(ref_data);

    auto slater = Slater<D>(alpha, coeff, pos);

    // Making ref value
    auto dist = math_utils::calc_distance<D>(slater.getPos(), ref);
    const auto ref_val = slater.getCoef() * std::exp(-slater.getAlpha() * std::abs(dist));

    // Making MRA
    const auto min_scale = -4;
    auto corner = std::array<int, D>{};
    auto boxes = std::array<int, D>{};
    corner.fill(-1);
    boxes.fill(2);
    auto world = BoundingBox<D>(min_scale, corner, boxes);
    const auto basis = InterpolatingBasis(7);
    auto MRA = MultiResolutionAnalysis<D>(world, basis, 25);

    // Making a function tree and projects the slater fcn onto it
    FunctionTree<D> f_tree(MRA);
    build_grid<D>(f_tree, slater);
    project(prec, f_tree, slater);
    f_tree.calcSquareNorm();

    WHEN("Slater function is projected") {
        THEN("The Slater function can be evaluated") { REQUIRE(slater.evalf(ref) == Catch::Approx(ref_val)); }
        THEN("The projected Slater function can be evaluated") { REQUIRE(f_tree.evalf(ref) == Catch::Approx(ref_val)); }
        THEN("the square norm matches the analytical value") { REQUIRE(f_tree.getSquareNorm() == Catch::Approx(slater.calcSquareNorm())); }
    }
}

    /*
SCENARIO("Slater_2D", "[slater_2d]") {
    const auto D = 2;
    const auto prec = 1.0e-4;
    // Making normalized gaussian centered at the origin
    const auto pos = Coord<D>{0.0, 0.0};
    const auto alpha = 2.2;
    const auto coeff = 3.3;   
    auto slater = Slater<D>(alpha, coeff, pos);

    // Making ref value
    const auto r_ref = std::array<double, D>{.1, 0.2};
    double dist = math_utils::calc_distance<D>(r_ref, slater.getPos());
    const auto ref_val = slater.getCoef() * std::exp(-slater.getAlpha() * std::abs(dist));

    // Making MRA
    const auto min_scale = -4;
    const auto corner = std::array<int, D>{-1,-1};
    const auto boxes = std::array<int, D>{2,2};
    auto world = BoundingBox<D>(min_scale, corner, boxes);
    const auto basis = InterpolatingBasis(7);
    auto MRA = MultiResolutionAnalysis<D>(world, basis, 25);

    // Making a function tree and projects the slater fcn onto it
    FunctionTree<D> f_tree(MRA);
    build_grid<D>(f_tree, slater);
    project(prec, f_tree, slater);
    f_tree.calcSquareNorm();

    WHEN("Slater function is projected") {
        THEN("The Slater function can be evaluated") { REQUIRE(slater.evalf(r_ref) == Catch::Approx(ref_val)); }

        THEN("The projected Slater function can be evaluated") { REQUIRE(f_tree.evalf(r_ref) == Catch::Approx(ref_val)); }

        THEN("the square norm matches") { REQUIRE(f_tree.getSquareNorm() == Catch::Approx(slater.calcSquareNorm())); }

    }
}
    */
    
} // namespace gaussians
