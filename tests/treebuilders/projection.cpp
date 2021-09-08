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

#include "catch.hpp"

#include "factory_functions.h"

#include "functions/function_utils.h"

#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"

using namespace mrcpp;

namespace projection {

template <int D> void testProjectFunction();
template <int D> void testProjectNarrowPeriodicGaussian();
template <int D> void testProjectWidePeriodicGaussian();

SCENARIO("Projecting Gaussian function", "[projection], [tree_builder], [trees]") {
    GIVEN("a Gaussian of unit charge in 1D") { testProjectFunction<1>(); }
    GIVEN("a Gaussian of unit charge in 2D") { testProjectFunction<2>(); }
    GIVEN("a Gaussian of unit charge in 3D") { testProjectFunction<3>(); }
    GIVEN("A periodic narrow Gaussian of unit charge in 1D") { testProjectNarrowPeriodicGaussian<1>(); }
    GIVEN("A periodic narrow Gaussian of unit charge in 2D") { testProjectNarrowPeriodicGaussian<2>(); }
    GIVEN("A periodic narrow Gaussian of unit charge in 3D") { testProjectNarrowPeriodicGaussian<3>(); }
    GIVEN("A periodic wide Gaussian of unit charge in 1D") { testProjectWidePeriodicGaussian<1>(); }
    GIVEN("A periodic wide Gaussian of unit charge in 2D") { testProjectWidePeriodicGaussian<2>(); }
    GIVEN("A periodic wide Gaussian of unit charge in 3D") { testProjectWidePeriodicGaussian<3>(); }
}

template <int D> void testProjectFunction() {
    GaussFunc<D> *func = nullptr;
    initialize(&func);
    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize(&mra);

    WHEN("the function is projected on the default grid") {
        FunctionTree<D> tree(*mra);
        project(-1.0, tree, *func);
        THEN("it integrates to approximately one") { REQUIRE(tree.integrate() == Approx(1.0).margin(1.0)); }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = tree.getSquareNorm();
            REQUIRE(dot(tree, tree) == Approx(norm));
        }
    }
    WHEN("the function is projected on an adapted grid") {
        FunctionTree<D> tree(*mra);
        build_grid(tree, *func);
        project(-1.0, tree, *func);
        THEN("it integrates to approximately one") { REQUIRE(tree.integrate() == Approx(1.0).epsilon(1.0e-3)); }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = tree.getSquareNorm();
            REQUIRE(dot(tree, tree) == Approx(norm));
        }
    }
    WHEN("the function is projected with guaranteed precision") {
        const double prec = 1.0e-4;
        FunctionTree<D> f_tree(*mra);
        project(prec, f_tree, *func);
        THEN("it integrates to approximately one") { REQUIRE(f_tree.integrate() == Approx(1.0).epsilon(1.0e-8)); }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = f_tree.getSquareNorm();
            REQUIRE(dot(f_tree, f_tree) == Approx(norm));
        }
        AND_WHEN("the function is projected on an identical grid") {
            FunctionTree<D> g_tree(*mra);
            copy_grid(g_tree, f_tree);

            project(-1.0, g_tree, *func);
            THEN("it integrates to the same value") {
                const double charge = f_tree.integrate();
                REQUIRE(g_tree.integrate() == Approx(charge));
            }
            THEN("the dot product with the original is equal to their squared norm") {
                const double norm = f_tree.getSquareNorm();
                REQUIRE(dot(g_tree, f_tree) == Approx(norm));
            }
        }
    }
    finalize(&mra);
    finalize(&func);
}

template <int D> void testProjectNarrowPeriodicGaussian() {
    const auto prec = 1.0e-4;

    auto period = std::array<double, D>{};
    period.fill(2.0); // Creating a world with period 2 in each direction
    auto periodic = true;

    GaussFunc<D> *func = nullptr;
    initialize<D>(&func);
    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize<D>(&mra, periodic, period);

    // Periodify Gaussian
    auto periodic_func = function_utils::periodify<D>(*func, period);

    FunctionTree<D> f_tree(*mra);
    build_grid<D>(f_tree, periodic_func);
    project<D>(prec, f_tree, periodic_func);
    REQUIRE(f_tree.integrate() == Approx(1.0));
}

template <int D> void testProjectWidePeriodicGaussian() {
    const auto prec = 1.0e-4;

    auto period = std::array<double, D>{};
    period.fill(2.0); // Creating a world with period 2 in each direction
    auto pos = Coord<D>{};
    pos.fill(1.0);

    auto alpha = 1.0;
    auto beta = std::pow(alpha / pi, static_cast<double>(D) / 2.0);
    auto func = GaussFunc<D>(alpha, beta, pos);
    auto periodic_func = function_utils::periodify<D>(func, period);

    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize<D>(&mra, true, period);

    FunctionTree<D> f_tree(*mra);
    build_grid<D>(f_tree, periodic_func);
    project<D>(prec, f_tree, periodic_func);
    REQUIRE(f_tree.integrate() == Approx(1.0));
}

TEST_CASE("Crop FunctionTree", "[function_tree], [crop]") {
    const int D = 3;
    const auto prec_1 = 1.0e-3;
    const auto prec_21 = 1.0e-5;
    const auto prec_22 = 1.0e-2;
    Printer::init(0);

    GaussFunc<D> *func = nullptr;
    initialize(&func);
    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize(&mra);

    FunctionTree<D> tree_1(*mra);
    project(prec_1, tree_1, *func);

    FunctionTree<D> tree_2(*mra);
    project(prec_21, tree_2, *func);

    const int nodes_1 = tree_1.getNNodes();
    const int nodes_21 = tree_2.getNNodes();
    REQUIRE(nodes_21 > nodes_1);

    tree_2.crop(prec_22);

    const int nodes_22 = tree_2.getNNodes();
    REQUIRE(nodes_22 < nodes_21);
}

} // namespace projection
