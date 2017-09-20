#include "catch.hpp"

#include "BandWidth.h"

using namespace Eigen;

namespace band_width {

TEST_CASE("BandWidth", "[band_width]") {
    const int depth = 10;
    BandWidth bw(depth);

    SECTION("Default BandWidth") {
        REQUIRE( (bw.getDepth() == depth) );
        for (int n = 0; n < 2*depth; n++) {
            REQUIRE( bw.isEmpty(n) );
            REQUIRE( (bw.getMaxWidth(n) == -1) );
        }
    }

    SECTION("Setting bandnd widths") {
        bw.setWidth(0, 0, 1);
        bw.setWidth(0, 1, 2);
        bw.setWidth(0, 2, 3);
        bw.setWidth(0, 3, 0);

        bw.setWidth(1, 0, 5);
        bw.setWidth(1, 1, 4);
        bw.setWidth(1, 2, 3);
        bw.setWidth(1, 3, 1);

        REQUIRE( bw.getDepth() == depth );
        REQUIRE_FALSE( bw.isEmpty(0) );
        REQUIRE_FALSE( bw.isEmpty(1) );
        REQUIRE( (bw.getMaxWidth(0) == 3) );
        REQUIRE( (bw.getMaxWidth(1) == 5) );
        for (int n = 2; n < 2*depth; n++) {
            REQUIRE( bw.isEmpty(n) );
            REQUIRE( (bw.getMaxWidth(n) == -1) );
        }
    }
}

} // namespace
