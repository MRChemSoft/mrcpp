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
        auto *cIdx = new NodeIndex<D>(*nIdx);
        testInitial<D>(cIdx);
        finalize(&cIdx);
    }

    SECTION("Child constructor") {
        int i = D;
        auto *cIdx = new NodeIndex<D>(*nIdx, i);
        REQUIRE( cIdx->getScale() == (nIdx->getScale() + 1) );
        finalize(&cIdx);
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

template<int D> void testCompare() {
    NodeIndex<D> *aIdx = 0;
    NodeIndex<D> *bIdx = 0;
    initialize<D>(&aIdx);
    initialize<D>(&bIdx);
    THEN("aIdx == bIdx") {
        REQUIRE( *aIdx == *bIdx );
        REQUIRE_FALSE( *aIdx != *bIdx );
    }
    WHEN("aIdx is given a deeper scale") {
        aIdx->setScale(2);
        THEN("aIdx != bIdx") {
            REQUIRE( *aIdx != *bIdx );
            REQUIRE_FALSE( *aIdx == *bIdx );
        }
    }
    WHEN("bIdx is given a coarser scale") {
        bIdx->setScale(-10);
        THEN("aIdx != bIdx") {
            REQUIRE( *aIdx != *bIdx );
            REQUIRE_FALSE( *aIdx == *bIdx );
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
            REQUIRE( *aIdx != *bIdx );
            REQUIRE_FALSE( *aIdx == *bIdx );
        }
    }
    finalize(&aIdx);
    finalize(&bIdx);
}

} // namespace
