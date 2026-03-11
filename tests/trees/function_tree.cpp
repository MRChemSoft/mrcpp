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

#include "treebuilders/multiply.h"

using namespace mrcpp;

namespace function_tree {

template <int D> void testZeroFunction();
template <int D> void testGeneratedNodes();

SCENARIO("Zero FunctionTree", "[function_tree_zero], [function_tree], [trees]") {
    GIVEN("a default function in 1D") { testZeroFunction<1>(); }
    GIVEN("a default function in 2D") { testZeroFunction<2>(); }
    GIVEN("a default function in 3D") { testZeroFunction<3>(); }
}

template <int D> void testZeroFunction() {
    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize(&mra);
    FunctionTree<D> tree(*mra);
    WHEN("the function is set to zero") {
        tree.setZero();
        THEN("its value in an arbitrary point is zero") {
            Coord<D> r;
            if (r.size() >= 1) r[0] = -0.20;
            if (r.size() >= 2) r[1] = 0.60;
            if (r.size() >= 3) r[2] = 0.76;
            REQUIRE(tree.evalf(r) == Catch::Approx(0.0));
        }
        THEN("its squared norm is zero") { REQUIRE(tree.getSquareNorm() == Catch::Approx(0.0)); }
        THEN("it integrates to zero") { REQUIRE(tree.integrate() == Catch::Approx(0.0)); }
        THEN("the dot product with itself is zero") { REQUIRE(dot(tree, tree) == Catch::Approx(0.0)); }
    }
    finalize(&mra);
}

SCENARIO("Generating FunctionTree nodes", "[function_tree_generating], [function_tree], [trees]") {
    GIVEN("a default function in 1D") { testGeneratedNodes<1>(); }
    GIVEN("a default function in 2D") { testGeneratedNodes<2>(); }
    GIVEN("a default function in 3D") { testGeneratedNodes<3>(); }
}

template <int D> void testGeneratedNodes() {
    unsigned int depth = 3;
    Coord<D> r;
    if (r.size() >= 1) r[0] = -0.3;
    if (r.size() >= 2) r[1] = 0.6;
    if (r.size() >= 3) r[2] = 1.9;

    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize(&mra);

    FunctionTree<D> tree(*mra);
    tree.setZero();

    THEN("there are no GenNodes") { REQUIRE(tree.getNGenNodes() == 0); }

    WHEN("a non-existing node is fetched") {
        MWNode<D> &node = tree.getNode(r, depth);
        THEN("there will be allocated GenNodes") {
            REQUIRE(tree.getNGenNodes() > 0);

            AND_WHEN("the GenNodes are deleted") {
                tree.deleteGenerated();
                THEN("there will be no GenNodes") { REQUIRE(tree.getNGenNodes() == 0); }
            }
        }
        (void)&node; // Clean up the fetched node
    }
    finalize(&mra);
}

} // namespace function_tree
