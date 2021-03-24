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
#include "trees/MWNode.h"

using namespace mrcpp;

namespace mw_tree {

template <int D> void testNodeFetchers();

TEST_CASE("MWTree: Fetching nodes", "[mw_tree_fetch], [mw_tree], [trees]") {
    SECTION("1D") { testNodeFetchers<1>(); }
    SECTION("2D") { testNodeFetchers<2>(); }
    SECTION("3D") { testNodeFetchers<3>(); }
}

template <int D> void testNodeFetchers() {
    Coord<D> r;
    if (D >= 1) r[0] = -0.3;
    if (D >= 2) r[1] = 0.6;
    if (D >= 3) r[2] = 1.9;

    const int cIdx = 1 << (D - 1);
    NodeIndex<D> *root = nullptr;
    initialize(&root);
    const NodeIndex<D> idx_0(*root);
    const NodeIndex<D> idx_1 = idx_0.child(cIdx);
    const NodeIndex<D> idx_2 = idx_1.child(cIdx);
    finalize(&root);

    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize(&mra);

    FunctionTree<D> tree(*mra);
    tree.setZero();

    const auto &const_tree = const_cast<const FunctionTree<D> &>(tree);

    // Fetch by NodeIndex
    SECTION("Find node by NodeIndex: existing node") {
        MWNode<D> *node = tree.findNode(idx_0);
        REQUIRE(node != 0);
        REQUIRE(node->getNodeIndex() == idx_0);
        REQUIRE(node->isAllocated());
        REQUIRE(node->hasCoefs());
        REQUIRE_FALSE(node->isGenNode());
    }
    SECTION("Find const node by NodeIndex: existing node") {
        const MWNode<D> *node = const_tree.findNode(idx_0);
        REQUIRE(node != 0);
        REQUIRE(node->getNodeIndex() == idx_0);
        REQUIRE(node->isAllocated());
        REQUIRE(node->hasCoefs());
        REQUIRE_FALSE(node->isGenNode());
    }
    SECTION("Find node by NodeIndex: non-existing node") {
        MWNode<D> *node = tree.findNode(idx_1);
        REQUIRE(node == 0);
    }
    SECTION("Find const node by NodeIndex: non-existing node") {
        const MWNode<D> *node = const_tree.findNode(idx_2);
        REQUIRE(node == 0);
    }
    SECTION("Get node by NodeIndex: existing node") {
        MWNode<D> &node = tree.getNode(idx_0);
        REQUIRE(node.getNodeIndex() == idx_0);
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE_FALSE(node.isGenNode());
    }
    SECTION("Get node by NodeIndex: non-existing node") {
        MWNode<D> &node = tree.getNode(idx_2);
        REQUIRE(node.getNodeIndex() == idx_2);
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.isGenNode());
    }
    SECTION("Get node or end node by NodeIndex: existing node") {
        MWNode<D> &node = tree.getNodeOrEndNode(idx_0);
        REQUIRE(node.getNodeIndex() == idx_0);
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE_FALSE(node.isGenNode());
    }
    SECTION("Get const node or end node by NodeIndex: existing node") {
        const MWNode<D> &node = const_tree.getNodeOrEndNode(idx_0);
        REQUIRE(node.getNodeIndex() == idx_0);
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE_FALSE(node.isGenNode());
    }
    SECTION("Get node or end node by NodeIndex: non-existing node") {
        MWNode<D> &node = tree.getNodeOrEndNode(idx_2);
        REQUIRE(node.getNodeIndex() != idx_2);
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.isEndNode());
        REQUIRE_FALSE(node.isGenNode());
    }
    SECTION("Get const node or end node by NodeIndex: non-existing node") {
        const MWNode<D> &node = const_tree.getNodeOrEndNode(idx_1);
        REQUIRE(node.getNodeIndex() != idx_2);
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.isEndNode());
        REQUIRE_FALSE(node.isGenNode());
    }

    // Fetch by coordinate
    SECTION("Get node by coord: existing node") {
        int depth = 0;
        MWNode<D> &node = tree.getNode(r, depth);
        REQUIRE(node.hasCoord(r));
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.getDepth() == depth);
    }
    SECTION("Get node by coord: non-existing node") {
        int depth = 3;
        MWNode<D> &node = tree.getNode(r, depth);
        REQUIRE(node.hasCoord(r));
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.getDepth() == depth);
    }
    SECTION("Get node or end node by coord: existing node") {
        int depth = 0;
        MWNode<D> &node = tree.getNodeOrEndNode(r, depth);
        REQUIRE(node.hasCoord(r));
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.getDepth() == depth);
    }
    SECTION("Get const node or end node by coord: existing node") {
        int depth = 0;
        const MWNode<D> &node = const_tree.getNodeOrEndNode(r, depth);
        REQUIRE(node.hasCoord(r));
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.getDepth() == depth);
    }
    SECTION("Get node or end node by coord: non-existing node") {
        int depth = 3;
        MWNode<D> &node = tree.getNodeOrEndNode(r, depth);
        REQUIRE(node.hasCoord(r));
        REQUIRE(node.isEndNode());
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.getDepth() != depth);
    }
    SECTION("Get const node or end node by coord: non-existing node") {
        int depth = 3;
        const MWNode<D> &node = const_tree.getNodeOrEndNode(r, depth);
        REQUIRE(node.hasCoord(r));
        REQUIRE(node.isEndNode());
        REQUIRE(node.isAllocated());
        REQUIRE(node.hasCoefs());
        REQUIRE(node.getDepth() != depth);
    }
    finalize(&mra);
}

} // namespace mw_tree
