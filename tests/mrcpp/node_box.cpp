#include "catch.hpp"

#include "factory_functions.h"

namespace node_box {

template<int D> void testConstructors();
template<int D> void testNodeFetchers();

TEST_CASE("NodeBox: Constructor", "[node_box_constructor], [node_box], [boxes]") {
    SECTION("1D") {
        testConstructors<1>();
    }
    SECTION("2D") {
        testConstructors<2>();
    }
    SECTION("3D") {
        testConstructors<3>();
    }
}

template<int D> void testConstructors() {
    int nb[D];
    int tot_boxes = 1;
    for (int d = 0; d < D; d++) {
        nb[d] = D+d;
        tot_boxes *= nb[d];
    }
    NodeIndex<D> *nIdx = 0;
    initialize(&nIdx);

    NodeBox<D> box(*nIdx, nb);
    finalize(&nIdx);

    SECTION("Constructor") {
        REQUIRE( (box.size() == tot_boxes) );
        REQUIRE( (box.getNOccupied() == 0) );
    }

    SECTION("Copy constructor") {
        NodeBox<D> box_copy(box);
        REQUIRE( (box_copy.size() == tot_boxes) );
        REQUIRE( (box_copy.getNOccupied() == 0) );
    }

    SECTION("Base class copy constructor") {
        const BoundingBox<D> &b_box = static_cast<const BoundingBox<D> &>(box);
        NodeBox<D> box_copy(b_box);
        REQUIRE( (box_copy == b_box) );
    }
}

TEST_CASE("NodeBox: Fetching nodes", "[node_box_fetch], [node_box], [boxes]") {
    SECTION("1D") {
        testNodeFetchers<1>();
    }
    SECTION("2D") {
        testNodeFetchers<2>();
    }
    SECTION("3D") {
        testNodeFetchers<3>();
    }
}

template<int D> void testNodeFetchers() {
    const double r[3] = {-0.3, 0.6, 1.9};

    int cIdx = 1 << (D - 1);
    NodeIndex<D> *root = 0;
    initialize(&root);
    const NodeIndex<D> idx_0(*root);
    const NodeIndex<D> idx_1(idx_0, cIdx);
    const NodeIndex<D> idx_2(idx_1, cIdx);
    finalize(&root);

    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);
    FunctionTree<D> tree(*mra);
    finalize(&mra);

    NodeBox<D> &node_box = tree.getRootBox();
    const NodeBox<D> &const_box = tree.getRootBox();

    // Fetch by NodeIndex
    SECTION("Find root node by NodeIndex") {
        MWNode<D> &node = node_box.getNode(idx_2);
        REQUIRE( (node.getDepth() == 0) );
        REQUIRE( node.isAncestor(idx_2) );
    }
    SECTION("Find const root node by NodeIndex") {
        const MWNode<D> &node = const_box.getNode(idx_2);
        REQUIRE( (node.getDepth() == 0) );
        REQUIRE( node.isAncestor(idx_2) );
    }

    // Fetch by coordinate
    SECTION("Get root node by coord: existing node") {
        MWNode<D> &node = node_box.getNode(r);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( (node.getDepth() == 0) );
    }
    SECTION("Get node by coord: non-existing node") {
        const MWNode<D> &node = const_box.getNode(r);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( (node.getDepth() == 0) );
    }
}

} // namespace
