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

#include "MRCPP/treebuilders/grid.h"
#include "MRCPP/treebuilders/multiply.h"
#include "MRCPP/treebuilders/project.h"

using namespace mrcpp;

namespace projection {

template <int D> void testProjectFunction();

SCENARIO("Projecting Gaussian function", "[projection], [tree_builder], [trees]") {
    GIVEN("a Gaussian of unit charge in 1D") { testProjectFunction<1>(); }
    GIVEN("a Gaussian of unit charge in 2D") { testProjectFunction<2>(); }
    GIVEN("a Gaussian of unit charge in 3D") { testProjectFunction<3>(); }
}

template <int D> void testProjectFunction() {
    GaussFunc<D> *func = 0;
    initialize(&func);
    MultiResolutionAnalysis<D> *mra = 0;
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

} // namespace projection
