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

#include "MRCPP/operators/ABGVOperator.h"
#include "MRCPP/operators/MWOperator.h"
#include "MRCPP/operators/PHOperator.h"
#include "MRCPP/treebuilders/add.h"
#include "MRCPP/treebuilders/apply.h"
#include "MRCPP/treebuilders/project.h"
#include "MRCPP/utils/math_utils.h"

using namespace mrcpp;

namespace derivative_operator {

template <int D> MultiResolutionAnalysis<D> *initializeMRA() {
    // Constructing world box
    int min_scale = -4;
    std::array<int, D> corner;
    std::array<int, D> boxes;
    std::array<double, D> scaling_factor;
    corner.fill(-1);
    boxes.fill(2);
    scaling_factor.fill(2.0);
    BoundingBox<D> world(min_scale, corner, boxes, scaling_factor);

    // Constructing scaling basis
    int order = 5;
    InterpolatingBasis basis(order);

    // Initializing MRA
    int max_depth = 20;
    return new MultiResolutionAnalysis<D>(world, basis, max_depth);
}

template <int D> MultiResolutionAnalysis<D> *initializePeriodicMRA() {
    // Constructing world box
    std::array<double, D> scaling_factor;
    scaling_factor.fill(2.0 * pi);
    auto periodic = true;
    BoundingBox<D> world(scaling_factor, periodic);

    // Constructing scaling basis
    int order = 5;
    InterpolatingBasis basis(order);

    // Initializing MRA
    int max_depth = 20;
    return new MultiResolutionAnalysis<D>(world, basis, max_depth);
}

template <int D> void testDifferentiationABGV(double a, double b) {
    MultiResolutionAnalysis<D> *mra = initializeMRA<D>();

    double prec = 1.0e-3;
    ABGVOperator<D> diff(*mra, a, b);

    Coord<D> r_0;
    for (auto &x : r_0) x = pi;

    auto f = [r_0](const Coord<D> &r) -> double {
        double R = math_utils::calc_distance<D>(r, r_0);
        return std::exp(-R * R);
    };

    auto df = [r_0](const Coord<D> &r) -> double {
        double R = math_utils::calc_distance<D>(r, r_0);
        return -2.0 * std::exp(-R * R) * (r[0] - r_0[0]);
    };

    FunctionTree<D> f_tree(*mra);
    project<D>(prec / 10, f_tree, f);

    FunctionTree<D> df_tree(*mra);
    project<D>(prec / 10, df_tree, df);

    FunctionTree<D> dg_tree(*mra);
    apply(dg_tree, diff, f_tree, 0);

    FunctionTree<D> err_tree(*mra);
    add(-1.0, err_tree, 1.0, df_tree, -1.0, dg_tree);

    double df_norm = std::sqrt(df_tree.getSquareNorm());
    double abs_err = std::sqrt(err_tree.getSquareNorm());
    double rel_err = abs_err / df_norm;

    REQUIRE(rel_err == Approx(0.0).margin(prec));

    delete mra;
}

template <int D> void testDifferentiationPH(int order) {
    MultiResolutionAnalysis<D> *mra = initializeMRA<D>();

    double prec = 1.0e-3;
    PHOperator<D> diff(*mra, order);

    Coord<D> r_0;
    for (auto &x : r_0) x = pi;

    auto f = [r_0](const Coord<D> &r) -> double {
        double R = math_utils::calc_distance<D>(r, r_0);
        return std::exp(-R * R);
    };
    auto df = [r_0, order](const Coord<D> &r) -> double {
        double R = math_utils::calc_distance<D>(r, r_0);
        return -(2 - order) * 2 * std::exp(-R * R) * (r[0] - pi) +
               (order - 1) * (-2 * std::exp(-R * R) + 4 * std::exp(-R * R) * (r[0] - r_0[0]) * (r[0] - r_0[0]));
        // 2-order = 1 and order-1 = 0 in the first order case
        // 2-order = 0 and order-1 = 1 in the second order case
    };

    FunctionTree<D> f_tree(*mra);
    project<D>(prec / 10, f_tree, f);

    FunctionTree<D> df_tree(*mra);
    project<D>(prec / 10, df_tree, df);

    FunctionTree<D> dg_tree(*mra);
    apply(dg_tree, diff, f_tree, 0);

    FunctionTree<D> err_tree(*mra);
    add(-1.0, err_tree, 1.0, df_tree, -1.0, dg_tree);

    double df_norm = std::sqrt(df_tree.getSquareNorm());
    double abs_err = std::sqrt(err_tree.getSquareNorm());
    double rel_err = abs_err / df_norm;

    REQUIRE(rel_err == Approx(0.0).margin(prec));

    delete mra;
}

template <int D> void testDifferentiationPeriodicABGV(double a, double b) {
    MultiResolutionAnalysis<D> *mra = initializePeriodicMRA<D>();

    double prec = 1.0e-3;
    ABGVOperator<D> diff(*mra, a, b);

    auto g_func = [](const mrcpp::Coord<D> &r) { return cos(r[0]) * cos(r[1]) * cos(r[2]); };
    auto dg_func = [](const mrcpp::Coord<D> &r) { return -sin(r[0]) * cos(r[1]) * cos(r[2]); };

    FunctionTree<D> g_tree(*mra);
    FunctionTree<D> dg_tree(*mra);

    project<D>(prec, g_tree, g_func);

    apply(dg_tree, diff, g_tree, 0);

    REQUIRE(dg_tree.evalf({0.0, 0.0, 0.0}) == Approx(dg_func({0.0, 0.0, 0.0})).margin(prec));
    REQUIRE(dg_tree.evalf({12.0, 0.0, 0.0}) == Approx(dg_func({12.0, 0.0, 0.0})).margin(prec));
    REQUIRE(dg_tree.evalf({20.0, 0.0, 0.0}) == Approx(dg_func({20.0, 0.0, 0.0})).margin(prec));
    REQUIRE(dg_tree.evalf({-5.0, 0.0, 0.0}) == Approx(dg_func({-5.0, 0.0, 0.0})).margin(prec));

    delete mra;
}

template <int D> void testDifferentiationPeriodicPH(int order) {
    MultiResolutionAnalysis<D> *mra = initializePeriodicMRA<D>();

    double prec = 1.0e-3;
    PHOperator<D> diff(*mra, order);

    auto g_func = [](const mrcpp::Coord<D> &r) { return cos(r[0]) * cos(r[1]) * cos(r[2]); };
    auto dg_func = [order](const mrcpp::Coord<D> &r) {
        if (order == 1) return -sin(r[0]) * cos(r[1]) * cos(r[2]); // First order
        return -cos(r[0]) * cos(r[1]) * cos(r[2]);                 // Second order
    };

    FunctionTree<D> g_tree(*mra);
    FunctionTree<D> dg_tree(*mra);

    project<D>(prec, g_tree, g_func);

    apply(dg_tree, diff, g_tree, 0);

    REQUIRE(dg_tree.evalf({0.0, 0.0, 0.0}) == Approx(dg_func({0.0, 0.0, 0.0})).margin(prec));
    REQUIRE(dg_tree.evalf({12.0, 0.0, 0.0}) == Approx(dg_func({12.0, 0.0, 0.0})).margin(prec));
    REQUIRE(dg_tree.evalf({20.0, 0.0, 0.0}) == Approx(dg_func({20.0, 0.0, 0.0})).margin(prec));
    REQUIRE(dg_tree.evalf({-5.0, 0.0, 0.0}) == Approx(dg_func({-5.0, 0.0, 0.0})).margin(prec));

    delete mra;
}

TEST_CASE("ABGV differentiantion central difference", "[derivative_operator], [central_difference]") {
    // 0.5,0.5 specifies central difference
    SECTION("1D derivative test") { testDifferentiationABGV<1>(0.5, 0.5); }
    SECTION("2D derivative test") { testDifferentiationABGV<2>(0.5, 0.5); }
    SECTION("3D derivative test") { testDifferentiationABGV<3>(0.5, 0.5); }
}

TEST_CASE("ABGV differentiantion center difference", "[derivative_operator], [center_difference]") {
    // 0,0 specifies center difference
    SECTION("1D derivative test") { testDifferentiationABGV<1>(0, 0); }
    SECTION("2D derivative test") { testDifferentiationABGV<2>(0, 0); }
    SECTION("3D derivative test") { testDifferentiationABGV<3>(0, 0); }
}

TEST_CASE("PH differentiantion first order", "[derivative_operator], [PH_first_order]") {
    SECTION("1D derivative test") { testDifferentiationPH<1>(1); }
    SECTION("2D derivative test") { testDifferentiationPH<2>(1); }
    SECTION("3D derivative test") { testDifferentiationPH<3>(1); }
}

TEST_CASE("PH differentiantion second order", "[derivative_operator], [PH_second_order]") {
    SECTION("1D second order derivative test") { testDifferentiationPH<1>(2); }
    SECTION("2D second order derivative test") { testDifferentiationPH<2>(2); }
    SECTION("3D second order derivative test") { testDifferentiationPH<3>(2); }
}

TEST_CASE("Periodic ABGV differentiantion central difference",
          "[periodic_derivative],[derivative_operator], [central_difference], [ABGV_periodic]") {
    // 0.5,0.5 specifies central difference
    SECTION("3D periodic derivative test") { testDifferentiationPeriodicABGV<3>(0.5, 0.5); }
}

TEST_CASE("Periodic PH differentiantion", "[periodic_derivative], [derivative_operator], [PH_periodic]") {
    // 0.5,0.5 specifies central difference
    SECTION("3D first order periodic derivative test") { testDifferentiationPeriodicPH<3>(1); }
    SECTION("3D first order periodic derivative test") { testDifferentiationPeriodicPH<3>(2); }
}

TEST_CASE("Gradient operator", "[derivative_operator], [gradient_operator]") {
    MultiResolutionAnalysis<3> *mra = initializeMRA<3>();

    double prec = 1.0e-3;
    ABGVOperator<3> diff(*mra, 0.0, 0.0);

    auto f = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return std::exp(-r2);
    };
    auto fx = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return -2.0 * r[0] * std::exp(-r2);
    };
    auto fy = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return -2.0 * r[1] * std::exp(-r2);
    };
    auto fz = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return -2.0 * r[2] * std::exp(-r2);
    };

    FunctionTree<3> f_tree(*mra);
    project<3>(prec, f_tree, f);

    FunctionTreeVector<3> grad_f = gradient(diff, f_tree);
    REQUIRE(grad_f.size() == 3);

    const Coord<3> r{1.1, 0.4, 0.2};
    REQUIRE(get_func(grad_f, 0).evalf(r) == Approx(fx(r)).epsilon(prec));
    REQUIRE(get_func(grad_f, 1).evalf(r) == Approx(fy(r)).epsilon(prec));
    REQUIRE(get_func(grad_f, 2).evalf(r) == Approx(fz(r)).epsilon(prec));
    clear(grad_f, true);

    delete mra;
}

TEST_CASE("Divergence operator", "[derivative_operator], [divergence_operator]") {
    MultiResolutionAnalysis<3> *mra = initializeMRA<3>();

    double prec = 1.0e-3;
    ABGVOperator<3> diff(*mra, 0.5, 0.5);

    auto f = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return std::exp(-r2);
    };
    auto fx = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return -2.0 * r[0] * std::exp(-r2);
    };
    auto fy = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return -2.0 * r[1] * std::exp(-r2);
    };
    auto fz = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return -2.0 * r[2] * std::exp(-r2);
    };

    FunctionTree<3> f_tree(*mra);
    project<3>(prec, f_tree, f);
    FunctionTreeVector<3> f_vec;
    f_vec.push_back(std::make_tuple(1.0, &f_tree));
    f_vec.push_back(std::make_tuple(2.0, &f_tree));
    f_vec.push_back(std::make_tuple(3.0, &f_tree));

    FunctionTree<3> div_f(*mra);
    divergence(div_f, diff, f_vec);

    for (int i = 0; i < 10; i++) {
        const Coord<3> r{-0.8 + 0.3 * i, 0.4 - 0.1 * i, 0.2 + 0.01 * i};
        const double ref = 1.0 * fx(r) + 2.0 * fy(r) + 3.0 * fz(r);
        REQUIRE(div_f.evalf(r) == Approx(ref).epsilon(prec));
    }

    delete mra;
}

} // namespace derivative_operator
