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
#include "trees/MWNode.h"

using namespace mrcpp;

namespace node_box {

template <int D> void testConstructors();
template <int D> void testNodeFetchers();

TEST_CASE("NodeBox: Constructor", "[node_box_constructor], [node_box], [boxes]") {
    SECTION("1D") { testConstructors<1>(); }
    SECTION("2D") { testConstructors<2>(); }
    SECTION("3D") { testConstructors<3>(); }
}

template <int D> void testConstructors() {
    std::array<int, D> nb;
    int tot_boxes = 1;
    for (int d = 0; d < D; d++) {
        nb[d] = D + d;
        tot_boxes *= nb[d];
    }
    NodeIndex<D> *nIdx = nullptr;
    initialize(&nIdx);

    NodeBox<D> box(*nIdx, nb);
    finalize(&nIdx);

    SECTION("Constructor") {
        REQUIRE(box.size() == tot_boxes);
        REQUIRE(box.getNOccupied() == 0);
    }

    SECTION("Copy constructor") {
        NodeBox<D> box_copy(box);
        REQUIRE(box_copy.size() == tot_boxes);
        REQUIRE(box_copy.getNOccupied() == 0);
    }

    SECTION("Base class copy constructor") {
        const auto &b_box = static_cast<const BoundingBox<D> &>(box);
        NodeBox<D> box_copy(b_box);
        REQUIRE(box_copy == b_box);
    }
}

TEST_CASE("NodeBox: Fetching nodes", "[node_box_fetch], [node_box], [boxes]") {
    SECTION("1D") { testNodeFetchers<1>(); }
    SECTION("2D") { testNodeFetchers<2>(); }
    SECTION("3D") { testNodeFetchers<3>(); }
}

template <int D> void testNodeFetchers() {
    Coord<D> r;
    if (D >= 1) r[0] = -0.3;
    if (D >= 2) r[1] = 0.6;
    if (D >= 3) r[2] = 1.9;

    int cIdx = 1 << (D - 1);
    NodeIndex<D> *root = nullptr;
    initialize(&root);
    const NodeIndex<D> idx_0(*root);
    const NodeIndex<D> idx_1(idx_0, cIdx);
    const NodeIndex<D> idx_2(idx_1, cIdx);
    finalize(&root);

    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize(&mra);
    FunctionTree<D> tree(*mra);
    finalize(&mra);

    NodeBox<D> &node_box = tree.getRootBox();
    const NodeBox<D> &const_box = tree.getRootBox();

    // Fetch by NodeIndex
    SECTION("Find root node by NodeIndex") {
        MWNode<D> &node = node_box.getNode(idx_2);
        REQUIRE(node.getDepth() == 0);
        REQUIRE(node.isAncestor(idx_2));
    }
    SECTION("Find const root node by NodeIndex") {
        const MWNode<D> &node = const_box.getNode(idx_2);
        REQUIRE(node.getDepth() == 0);
        REQUIRE(node.isAncestor(idx_2));
    }

    // Fetch by coordinate
    SECTION("Get root node by coord: existing node") {
        MWNode<D> &node = node_box.getNode(r);
        REQUIRE(node.hasCoord(r));
        REQUIRE(node.getDepth() == 0);
    }
    SECTION("Get node by coord: non-existing node") {
        const MWNode<D> &node = const_box.getNode(r);
        REQUIRE(node.hasCoord(r));
        REQUIRE(node.getDepth() == 0);
    }
}

} // namespace node_box
