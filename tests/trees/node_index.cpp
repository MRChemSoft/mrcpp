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

using namespace mrcpp;

namespace node_index {

template <int D> void testConstructors();
template <int D> void testRelations();
template <int D> void testCompare();

TEST_CASE("NodeIndex constructors", "[node_index_constructor], [node_index]") {
    SECTION("1D") { testConstructors<1>(); }
    SECTION("2D") { testConstructors<2>(); }
    SECTION("3D") { testConstructors<3>(); }
}

TEST_CASE("NodeIndex family relations", "[node_index_relations], [node_index]") {
    SECTION("1D") { testRelations<1>(); }
    SECTION("2D") { testRelations<2>(); }
    SECTION("3D") { testRelations<3>(); }
}

SCENARIO("Node indices can be compared", "[node_index_compare], [node_index]") {
    GIVEN("Identical node indices in 1D") { testCompare<1>(); }
    GIVEN("Identical node indices in 2D") { testCompare<2>(); }
    GIVEN("Identical node indices in 3D") { testCompare<3>(); }
}

template <int D> void testConstructors() {
    NodeIndex<D> *nIdx = nullptr;
    initialize<D>(&nIdx);

    SECTION("Constructor") { testInitial<D>(nIdx); }

    SECTION("Copy constructor") {
        auto *cIdx = new NodeIndex<D>(*nIdx);
        testInitial<D>(cIdx);
        finalize(&cIdx);
    }

    SECTION("Child constructor") {
        int i = D;
        auto cIdx = nIdx->child(i);
        REQUIRE(cIdx.getScale() == (nIdx->getScale() + 1));
    }

    SECTION("Parent constructor") {
        auto cIdx = nIdx->parent();
        REQUIRE(cIdx.getScale() == (nIdx->getScale() - 1));
    }
    SECTION("Default constructor") {
        auto *cIdx = new NodeIndex<D>();
        SECTION("Assignment operator") {
            *cIdx = *nIdx;
            testInitial<D>(cIdx);
        }
        finalize(&cIdx);
    }
    finalize(&nIdx);
}

template <int D> void testRelations() {
    for (int n = -2; n < 3; n++) {
        for (int l = -2; l < 5; l++) {
            NodeIndex<D> parent(n, {l});
            auto grand = parent.parent();
            auto child_0 = parent.child(D);
            auto child_00 = child_0.child(D);
            auto child_01 = child_0.child(D - 1);
            auto child_1 = parent.child(D - 1);
            auto child_10 = child_1.child(D);
            auto child_11 = child_1.child(D - 1);
            REQUIRE(related<D>(grand, grand));
            REQUIRE(related<D>(parent, grand));
            REQUIRE(related<D>(grand, child_0));
            REQUIRE(related<D>(child_1, grand));
            REQUIRE(related<D>(child_00, grand));
            REQUIRE(related<D>(grand, child_01));
            REQUIRE(related<D>(grand, child_10));
            REQUIRE(related<D>(child_11, grand));

            REQUIRE(related<D>(parent, grand));
            REQUIRE(related<D>(parent, parent));
            REQUIRE(related<D>(child_0, parent));
            REQUIRE(related<D>(parent, child_1));
            REQUIRE(related<D>(child_00, parent));
            REQUIRE(related<D>(parent, child_01));
            REQUIRE(related<D>(child_10, parent));
            REQUIRE(related<D>(parent, child_11));

            REQUIRE(related<D>(grand, child_0));
            REQUIRE(related<D>(child_0, parent));
            REQUIRE(!related<D>(child_1, child_0));
            REQUIRE(related<D>(child_0, child_01));
            REQUIRE(related<D>(child_00, child_0));
            REQUIRE(!related<D>(child_10, child_0));
            REQUIRE(!related<D>(child_0, child_11));

            REQUIRE(related<D>(child_1, grand));
            REQUIRE(related<D>(parent, child_1));
            REQUIRE(!related<D>(child_0, child_1));
            REQUIRE(related<D>(child_11, child_1));
            REQUIRE(!related<D>(child_1, child_00));
            REQUIRE(related<D>(child_0, child_00));
            REQUIRE(!related<D>(child_01, child_1));

            REQUIRE(siblings<D>(child_1, child_0));
            REQUIRE(siblings<D>(child_01, child_00));
            REQUIRE(!siblings<D>(child_01, child_11));
            REQUIRE(!siblings<D>(child_11, child_1));
            REQUIRE(!siblings<D>(grand, child_1));
            REQUIRE(!siblings<D>(grand, parent));
        }
    }
}

template <int D> void testCompare() {
    NodeIndex<D> *aIdx = nullptr;
    NodeIndex<D> *bIdx = nullptr;
    initialize<D>(&aIdx);
    initialize<D>(&bIdx);
    THEN("aIdx == bIdx") {
        REQUIRE(*aIdx == *bIdx);
        REQUIRE_FALSE(*aIdx != *bIdx);
    }
    WHEN("aIdx is given a deeper scale") {
        aIdx->setScale(2);
        THEN("aIdx != bIdx") {
            REQUIRE(*aIdx != *bIdx);
            REQUIRE_FALSE(*aIdx == *bIdx);
        }
    }
    WHEN("bIdx is given a coarser scale") {
        bIdx->setScale(-10);
        THEN("aIdx != bIdx") {
            REQUIRE(*aIdx != *bIdx);
            REQUIRE_FALSE(*aIdx == *bIdx);
        }
    }
    WHEN("aIdx is given a different translation") {
        (*aIdx)[D - 1]++;
        THEN("aIdx != bIdx") {
            REQUIRE(*aIdx != *bIdx);
            REQUIRE_FALSE(*aIdx == *bIdx);
        }
    }
    finalize(&aIdx);
    finalize(&bIdx);
}

} // namespace node_index
