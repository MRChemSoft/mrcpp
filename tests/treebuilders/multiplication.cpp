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

#include "functions/GaussPoly.h"
#include "treebuilders/WaveletAdaptor.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"

using namespace mrcpp;

namespace multiplication {

template <int D> void testMultiplication();
template <int D> void testSquare();

SCENARIO("Multiplying MW trees", "[multiplication], [tree_builder]") {
    GIVEN("Two MW functions in 1D") { testMultiplication<1>(); }
    GIVEN("Two MW functions in 2D") { testMultiplication<2>(); }
    GIVEN("Two MW functions in 3D") { testMultiplication<3>(); }
}

template <int D> void testMultiplication() {
    const double prec = 1.0e-4;

    double alpha = 1.0;
    double beta_a = 110.0;
    double beta_b = 50.0;

    double pos_c_a[3] = {-0.25, 0.35, 1.05};
    auto pos_a = details::convert_to_std_array<double, D>(pos_c_a);
    double pos_c_b[3] = {-0.20, 0.50, 1.05};
    auto pos_b = details::convert_to_std_array<double, D>(pos_c_b);

    GaussFunc<D> a_func(beta_a, alpha, pos_a);
    GaussFunc<D> b_func(beta_b, alpha, pos_b);
    GaussPoly<D> ref_func = a_func * b_func;

    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    // Initialize trees
    FunctionTree<D> a_tree(*mra);
    FunctionTree<D> b_tree(*mra);
    FunctionTree<D> ref_tree(*mra);

    // Build empty grids
    build_grid(a_tree, a_func);
    build_grid(b_tree, b_func);
    build_grid(ref_tree, ref_func);

    // Project functions
    project(prec, a_tree, a_func);
    project(prec, b_tree, b_func);
    project(prec, ref_tree, ref_func);

    const double ref_int = ref_tree.integrate();
    const double ref_norm = ref_tree.getSquareNorm();

    FunctionTreeVector<D> prod_vec;
    WHEN("the functions are multiplied") {
        FunctionTree<D> c_tree(*mra);
        prod_vec.push_back(std::make_tuple(1.0, &a_tree));
        prod_vec.push_back(std::make_tuple(1.0, &b_tree));
        multiply(prec, c_tree, prod_vec);
        prod_vec.clear();

        THEN("the MW product equals the analytic product") {
            double c_int = c_tree.integrate();
            double c_dot = dot(c_tree, ref_tree);
            double c_norm = c_tree.getSquareNorm();
            REQUIRE(c_int == Approx(ref_int));
            REQUIRE(c_dot == Approx(ref_norm));
            REQUIRE(c_norm == Approx(ref_norm));
        }
    }
    WHEN("the functions are multiplied in-place") {
        a_tree.multiply(1.0, b_tree);

        THEN("the first function equals the analytic product") {
            double a_int = a_tree.integrate();
            double a_dot = dot(a_tree, ref_tree);
            double a_norm = a_tree.getSquareNorm();
            REQUIRE(a_int == Approx(ref_int));
            REQUIRE(a_dot == Approx(ref_norm));
            REQUIRE(a_norm == Approx(ref_norm));
        }
    }
    finalize(&mra);
}

SCENARIO("Squaring MW trees", "[square], [tree_builder]") {
    GIVEN("A MW function in 1D") { testSquare<1>(); }
    GIVEN("A MW function in 2D") { testSquare<2>(); }
    GIVEN("A MW function in 3D") { testSquare<3>(); }
}

template <int D> void testSquare() {
    const double prec = 1.0e-4;

    double alpha = 1.0;
    double beta = 50.0;
    double pos_c[3] = {-0.25, 0.35, 1.05};
    auto pos = details::convert_to_std_array<double, D>(pos_c);

    GaussFunc<D> f_func(beta, alpha, pos);
    GaussPoly<D> ref_func = f_func * f_func;

    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    // Initialize trees
    FunctionTree<D> f_tree(*mra);
    FunctionTree<D> ref_tree(*mra);

    // Build empty grids
    build_grid(f_tree, f_func);
    build_grid(ref_tree, ref_func);

    // Project functions
    project(prec, f_tree, f_func);
    project(prec, ref_tree, ref_func);

    const double ref_int = ref_tree.integrate();
    const double ref_norm = ref_tree.getSquareNorm();

    WHEN("the function is squared out-of-place") {
        FunctionTree<D> ff_tree(*mra);
        square(prec, ff_tree, f_tree);

        THEN("the MW square equals the analytic square") {
            double ff_int = ff_tree.integrate();
            double ff_dot = dot(ff_tree, ref_tree);
            double ff_norm = ff_tree.getSquareNorm();
            REQUIRE(ff_int == Approx(ref_int));
            REQUIRE(ff_dot == Approx(ref_norm));
            REQUIRE(ff_norm == Approx(ref_norm));
        }
    }
    WHEN("the function is squared in-place") {
        f_tree.square();

        THEN("the MW square equals the analytic square") {
            double f_int = f_tree.integrate();
            double f_dot = dot(f_tree, ref_tree);
            double f_norm = f_tree.getSquareNorm();
            REQUIRE(f_int == Approx(ref_int));
            REQUIRE(f_dot == Approx(ref_norm));
            REQUIRE(f_norm == Approx(ref_norm));
        }
    }
    WHEN("the function is raised to the second power out-of-place") {
        FunctionTree<D> ff_tree(*mra);
        power(prec, ff_tree, f_tree, 2.0);

        THEN("the MW power equals the analytic power") {
            double ff_int = ff_tree.integrate();
            double ff_dot = dot(ff_tree, ref_tree);
            double ff_norm = ff_tree.getSquareNorm();
            REQUIRE(ff_int == Approx(ref_int));
            REQUIRE(ff_dot == Approx(ref_norm));
            REQUIRE(ff_norm == Approx(ref_norm));
        }
    }
    WHEN("the function is raised to second power in-place") {
        f_tree.power(2.0);

        THEN("the MW power(2.0) equals the analytic square") {
            double f_int = f_tree.integrate();
            double f_dot = dot(f_tree, ref_tree);
            double f_norm = f_tree.getSquareNorm();
            REQUIRE(f_int == Approx(ref_int));
            REQUIRE(f_dot == Approx(ref_norm));
            REQUIRE(f_norm == Approx(ref_norm));
        }
    }
    finalize(&mra);
}

TEST_CASE("Dot product FunctionTreeVectors", "[multiplication], [tree_vector_dot]") {
    MultiResolutionAnalysis<3> *mra = nullptr;
    initialize(&mra);

    double prec = 1.0e-4;

    auto fx = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return r[1] * r[2] * std::exp(-1.0 * r2);
    };
    auto fy = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return r[0] * r[2] * std::exp(-1.5 * r2);
    };
    auto fz = [](const Coord<3> &r) -> double {
        double r2 = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return r[0] * r[1] * std::exp(-2.0 * r2);
    };

    FunctionTree<3> fx_tree(*mra);
    FunctionTree<3> fy_tree(*mra);
    FunctionTree<3> fz_tree(*mra);

    project<3>(prec, fx_tree, fx);
    project<3>(prec, fy_tree, fy);
    project<3>(prec, fz_tree, fz);

    FunctionTreeVector<3> vec_a;
    vec_a.push_back(std::make_tuple(1.0, &fx_tree));
    vec_a.push_back(std::make_tuple(2.0, &fy_tree));
    vec_a.push_back(std::make_tuple(3.0, &fz_tree));

    FunctionTreeVector<3> vec_b;
    vec_b.push_back(std::make_tuple(1.0, &fz_tree));
    vec_b.push_back(std::make_tuple(2.0, &fy_tree));
    vec_b.push_back(std::make_tuple(3.0, &fx_tree));

    FunctionTree<3> dot_ab(*mra);
    build_grid(dot_ab, vec_a);
    build_grid(dot_ab, vec_b);
    dot(0.1 * prec, dot_ab, vec_a, vec_b);

    for (int i = 0; i < 10; i++) {
        const Coord<3> r = {-0.4 + 0.01 * i, 0.9 - 0.05 * i, 0.7 + 0.1 * i};
        const double ref = 1.0 * fx(r) * fz(r) + 4.0 * fy(r) * fy(r) + 9.0 * fz(r) * fx(r);
        REQUIRE(dot_ab.evalf(r) == Approx(ref).epsilon(prec));
    }

    finalize(&mra);
}

} // namespace multiplication
