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

using namespace mrcpp;

namespace bounding_box {

template <int D> void testConstructors();
template <int D> void testCompare();
template <int D> void testFetch();

TEST_CASE("BoundingBox constructors", "[bounding_box_constructor], [bounding_box], [boxes]") {
    SECTION("1D") { testConstructors<1>(); }
    SECTION("2D") { testConstructors<2>(); }
    SECTION("3D") { testConstructors<3>(); }
}

TEST_CASE("Bounding box comparison", "[bounding_box_compare], [bounding_box], [boxes]") {
    SECTION("1D") { testCompare<1>(); }
    SECTION("2D") { testCompare<2>(); }
    SECTION("3D") { testCompare<3>(); }
}

TEST_CASE("Bounding box fetching", "[bounding_box_fetch], [bounding_box], [boxes]") {
    SECTION("1D") { testFetch<1>(); }
    SECTION("2D") { testFetch<2>(); }
    SECTION("3D") { testFetch<3>(); }
}

template <int D> void testConstructors() {
    BoundingBox<D> *box = nullptr;
    initialize<D>(&box);

    SECTION("Constructor") { testInitial<D>(box); }

    SECTION("Copy constructor") {
        auto *box_copy = new BoundingBox<D>(*box);
        testInitial<D>(box_copy);
        finalize(&box_copy);
    }

    SECTION("Default constructor") {
        auto *box_copy = new BoundingBox<D>();
        SECTION("Assignment operator") {
            *box_copy = *box;
            testInitial<D>(box);
        }
        finalize(&box_copy);
    }
    finalize(&box);
}

template <int D> void testCompare() {
    SECTION("Identical boxes") {
        BoundingBox<D> *aBox = nullptr;
        BoundingBox<D> *bBox = nullptr;
        initialize<D>(&aBox);
        initialize<D>(&bBox);
        REQUIRE(*aBox == *bBox);
        REQUIRE_FALSE(*aBox != *bBox);
        finalize(&aBox);
        finalize(&bBox);
    }
    SECTION("Different boxes") {
        auto *aBox = new BoundingBox<D>();
        BoundingBox<D> *bBox = nullptr;
        initialize<D>(&bBox);
        REQUIRE(*aBox != *bBox);
        REQUIRE_FALSE(*aBox == *bBox);
        finalize(&aBox);
        finalize(&bBox);
    }
}

template <int D> void testFetch() {
    BoundingBox<D> *box = nullptr;
    initialize(&box);
    SECTION("Fetch by coord") {
        Coord<D> r;
        SECTION("Within bounds") {
            for (int d = 0; d < D; d++) { r[d] = box->getUpperBound(d) - 1.0e-15; }
            const int last = box->size() - 1;
            REQUIRE(box->getBoxIndex(r) == last);
        }
        SECTION("Out of bounds") {
            for (int d = 0; d < D; d++) { r[d] = -1.0; }
            REQUIRE(box->getBoxIndex(r) == -1);
        }
    }
    SECTION("Fetch by index") {
        SECTION("Within bounds") {
            NodeIndex<D> *idx = nullptr;
            initialize(&idx);
            REQUIRE(box->getBoxIndex(*idx) == 0);
            finalize(&idx);
        }
        SECTION("Out of bounds") {
            int n = -10;
            NodeIndex<D> idx(n);
            REQUIRE(box->getBoxIndex(idx) < 0);
        }
    }
    finalize(&box);
}

} // namespace bounding_box
