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

#include "treebuilders/project.h"
#include "treebuilders/grid.h"

using namespace mrcpp;

namespace grid_clear {

template<int D> void testClearGrid();

SCENARIO("Projected trees can be cleared and reused", "[clear_grid], [tree_builder], [trees]") {
    GIVEN("A projected FunctionTree in 1D") {
        testClearGrid<1>();
    }
    GIVEN("A projected FunctionTree in 2D") {
        testClearGrid<2>();
    }
    GIVEN("A projected FunctionTree in 3D") {
        testClearGrid<3>();
    }
}

template<int D> void testClearGrid() {
    GaussFunc<D> *func = 0;
    initialize(&func);
    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    const double prec = 1.0e-4;
    FunctionTree<D> tree(*mra);
    build_grid(tree, *func);
    project(prec, tree, *func);

    const int refDepth = tree.getDepth();
    const int refNodes = tree.getNNodes();
    const double refInt = tree.integrate();
    const double refNorm = tree.getSquareNorm();

    WHEN("the tree is cleared") {
        tree.clear();
        THEN("it represents an undefined function on the root grid") {
            REQUIRE( tree.getDepth() == 1 );
            REQUIRE( tree.getNNodes() == tree.getRootBox().size() );
            REQUIRE( tree.integrate() == Approx(0.0) );
            REQUIRE( tree.getSquareNorm() == Approx(-1.0) );
            AND_WHEN("the function is re-projected") {
                build_grid(tree, *func);
                project(prec, tree, *func);
                THEN("the representation is the same as before it was cleared") {
                    REQUIRE( tree.getDepth() == refDepth );
                    REQUIRE( tree.getNNodes() == refNodes );
                    REQUIRE( tree.integrate() == Approx(refInt) );
                    REQUIRE( tree.getSquareNorm() == Approx(refNorm) );
                }
            }
        }
    }
    WHEN("the grid is cleared") {
        clear_grid(tree);
        THEN("it represents an undefined function on the same grid") {
            REQUIRE( tree.getDepth() == refDepth );
            REQUIRE( tree.getNNodes() == refNodes );
            REQUIRE( tree.integrate() == Approx(0.0) );
            REQUIRE( tree.getSquareNorm() == Approx(-1.0) );
            AND_WHEN("the function is re-projected on the same grid") {
                project(-1.0, tree, *func);
                THEN("the representation is the same as before") {
                    REQUIRE( tree.getDepth() == refDepth );
                    REQUIRE( tree.getNNodes() == refNodes );
                    REQUIRE( tree.integrate() == Approx(refInt) );
                    REQUIRE( tree.getSquareNorm() == Approx(refNorm) );
                }
            }
        }
    }
    WHEN("the grid is refined adaptively") {
        const double new_prec = 1.0e-5;
        refine_grid(tree, new_prec);
        THEN("it represents the same function on a larger grid") {
            REQUIRE( tree.getDepth() >= refDepth );
            REQUIRE( tree.getNNodes() > refNodes );
            REQUIRE( tree.integrate() == Approx(refInt).epsilon(1.0e-8) );
            REQUIRE( tree.getSquareNorm() == Approx(refNorm).epsilon(1.0e-8) );
        }
    }
    finalize(&mra);
    finalize(&func);
}

} // namespace clear_grid
