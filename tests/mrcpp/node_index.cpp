#include "catch.hpp"

#include "factory_functions.h"

namespace node_index {

template<int D> void testConstructors();
template<int D> void testCompare();

TEST_CASE("NodeIndex constructors", "[node_index_constructor], [node_index]") {
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

SCENARIO("Node indices can be compared", "[node_index_compare], [node_index]") {
    GIVEN("Identical node indices in 1D") {
        testCompare<1>();
    }
    GIVEN("Identical node indices in 2D") {
        testCompare<2>();
    }
    GIVEN("Identical node indices in 3D") {
        testCompare<3>();
    }
}

template<int D> void testConstructors() {
    NodeIndex<D> *nIdx = 0;
    initialize<D>(&nIdx);

    SECTION("Constructor") {
        testInitial<D>(nIdx);
    }

    SECTION("Copy constructor") {
        NodeIndex<D> *cIdx = new NodeIndex<D>(*nIdx);
        testInitial<D>(cIdx);
        finalize(&cIdx);
    }

    SECTION("Child constructor") {
        int i = D;
        NodeIndex<D> *cIdx = new NodeIndex<D>(*nIdx, i);
        REQUIRE( (cIdx->getScale() == (nIdx->getScale() + 1)) );
        finalize(&cIdx);
    }

    SECTION("Default constructor") {
        NodeIndex<D> *cIdx = new NodeIndex<D>();
        SECTION("Assignment operator") {
            *cIdx = *nIdx;
            testInitial<D>(cIdx);
        }
        finalize(&cIdx);
    }
    finalize(&nIdx);
}

template<int D> void testCompare() {
    NodeIndex<D> *aIdx = 0;
    NodeIndex<D> *bIdx = 0;
    initialize<D>(&aIdx);
    initialize<D>(&bIdx);
    THEN("aIdx == bIdx") {
        REQUIRE( (*aIdx == *bIdx) );
        REQUIRE_FALSE( (*aIdx != *bIdx) );
    }
    WHEN("aIdx is given a deeper scale") {
        aIdx->setScale(2);
        THEN("aIdx != bIdx") {
            REQUIRE( (*aIdx != *bIdx) );
            REQUIRE_FALSE( (*aIdx == *bIdx) );
        }
    }
    WHEN("bIdx is given a coarser scale") {
        bIdx->setScale(-10);
        THEN("aIdx != bIdx") {
            REQUIRE( (*aIdx != *bIdx) );
            REQUIRE_FALSE( (*aIdx == *bIdx) );
        }
    }
    WHEN("aIdx is given a different translation") {
        const int *l1 = aIdx->getTranslation();
        int l2[D];
        for (int d = 0; d < D; d++) {
            l2[d] = l1[d];
        }
        l2[D-1]++;
        aIdx->setTranslation(l2);
        THEN("aIdx != bIdx") {
            REQUIRE( (*aIdx != *bIdx) );
            REQUIRE_FALSE( (*aIdx == *bIdx) );
        }
    }
    finalize(&aIdx);
    finalize(&bIdx);
}

} // namespace
