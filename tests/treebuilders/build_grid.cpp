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

using namespace mrcpp;

namespace grid_build {

template <int D> void testBuildGrid();

SCENARIO("Building empty grids", "[build_grid], [tree_builder], [trees]") {
    GIVEN("an analytic function in 1D") { testBuildGrid<1>(); }
    GIVEN("an analytic function in 2D") { testBuildGrid<2>(); }
    GIVEN("an analytic function in 3D") { testBuildGrid<3>(); }
}

template <int D> void testBuildGrid() {
    GaussFunc<D> *f_func = nullptr;
    initialize(&f_func);

    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize(&mra);

    WHEN("the GridGenerator is given no argument") {
        FunctionTree<D> f_tree(*mra);

        THEN("we get an empty tree of root nodes") {
            REQUIRE(f_tree.getSquareNorm() == Approx(-1.0));
            REQUIRE(f_tree.getDepth() == 1);
            REQUIRE(f_tree.getNNodes() == f_tree.getNEndNodes());
            REQUIRE(f_tree.getNGenNodes() == 0);

            AND_WHEN("the GridGenerator is given the analytic function") {
                build_grid(f_tree, *f_func, 2);

                THEN("the empty tree gets adapted") {
                    REQUIRE(f_tree.getSquareNorm() == Approx(-1.0));
                    REQUIRE(f_tree.getDepth() == 3);
                    REQUIRE(f_tree.getNNodes() > f_tree.getNEndNodes());
                    REQUIRE(f_tree.getNGenNodes() == 0);
                }
            }
        }
    }

    WHEN("the GridGenerator is given the analytic function") {
        FunctionTree<D> f_tree(*mra);
        build_grid(f_tree, *f_func, 2);

        THEN("we get an empty adapted tree structure") {
            REQUIRE(f_tree.getSquareNorm() == Approx(-1.0));
            REQUIRE(f_tree.getDepth() == 3);
            REQUIRE(f_tree.getNNodes() > f_tree.getNEndNodes());
            REQUIRE(f_tree.getNGenNodes() == 0);

            AND_WHEN("the empty tree is passed to the GridGenerator") {
                FunctionTree<D> g_tree(*mra);
                build_grid(g_tree, f_tree);

                THEN("we get an identical empty grid") {
                    REQUIRE(g_tree.getSquareNorm() == Approx(-1.0));
                    REQUIRE(g_tree.getDepth() == f_tree.getDepth());
                    REQUIRE(g_tree.getNNodes() == f_tree.getNNodes());
                    REQUIRE(g_tree.getNEndNodes() == f_tree.getNEndNodes());
                    REQUIRE(g_tree.getNGenNodes() == 0);
                }
            }
        }
    }
    finalize(&mra);
    finalize(&f_func);
}

} // namespace grid_build
