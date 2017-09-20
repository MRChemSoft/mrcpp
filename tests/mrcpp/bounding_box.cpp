#include "catch.hpp"

#include "factory_functions.h"

namespace bounding_box {

template<int D> void testConstructors();
template<int D> void testCompare();
template<int D> void testFetch();

TEST_CASE("BoundingBox constructors", "[bounding_box_constructor], [bounding_box], [boxes]") {
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

TEST_CASE("Bounding box comparison", "[bounding_box_compare], [bounding_box], [boxes]") {
    SECTION("1D") {
        testCompare<1>();
    }
    SECTION("2D") {
        testCompare<2>();
    }
    SECTION("3D") {
        testCompare<3>();
    }
}

TEST_CASE("Bounding box fetching", "[bounding_box_fetch], [bounding_box], [boxes]") {
    SECTION("1D") {
        testFetch<1>();
    }
    SECTION("2D") {
        testFetch<2>();
    }
    SECTION("3D") {
        testFetch<3>();
    }
}

template<int D> void testConstructors() {
    BoundingBox<D> *box = 0;
    initialize<D>(&box);

    SECTION("Constructor") {
        testInitial<D>(box);
    }

    SECTION("Copy constructor") {
        BoundingBox<D> *box_copy = new BoundingBox<D>(*box);
        testInitial<D>(box_copy);
        finalize(&box_copy);
    }

    SECTION("Default constructor") {
        BoundingBox<D> *box_copy = new BoundingBox<D>();
        SECTION("Assignment operator") {
            *box_copy = *box;
            testInitial<D>(box);
        }
        finalize(&box_copy);
    }
    finalize(&box);
}

template<int D> void testCompare() {
    SECTION("Identical boxes") {
        BoundingBox<D> *aBox = 0;
        BoundingBox<D> *bBox = 0;
        initialize<D>(&aBox);
        initialize<D>(&bBox);
        REQUIRE( (*aBox == *bBox) );
        REQUIRE_FALSE( (*aBox != *bBox) );
        finalize(&aBox);
        finalize(&bBox);
    }
    SECTION("Different boxes") {
        BoundingBox<D> *aBox = new BoundingBox<D>();
        BoundingBox<D> *bBox = 0;
        initialize<D>(&bBox);
        REQUIRE( (*aBox != *bBox) );
        REQUIRE_FALSE( (*aBox == *bBox) );
        finalize(&aBox);
        finalize(&bBox);
    }
}

template<int D> void testFetch() {
    BoundingBox<D> *box = 0;
    initialize(&box);
    SECTION("Fetch by coord") {
        double r[D];
        SECTION("Within bounds") {
            for (int d = 0; d < D; d++) {
                r[d] = box->getUpperBound(d) - 1.0e-15;
            }
            const int last = box->size() - 1;
            REQUIRE( (box->getBoxIndex(r) == last) );
        }
        SECTION("Out of bounds") {
            for (int d = 0; d < D; d++) {
                r[d] = -1.0;
            }
            REQUIRE( (box->getBoxIndex(r) == -1) );
        }
    }
    SECTION("Fetch by index") {
        SECTION("Within bounds") {
            NodeIndex<D> *idx = 0;
            initialize(&idx);
            REQUIRE( (box->getBoxIndex(*idx) == 0) );
            finalize(&idx);
        }
        SECTION("Out of bounds") {
            int n = -10;
            NodeIndex<D> idx(n);
            REQUIRE( (box->getBoxIndex(idx) < 0) );
        }
    }
    finalize(&box);
}

}
