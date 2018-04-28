#include "catch.hpp"

#include "core/FilterCache.h"
#include "constants.h"

using namespace Eigen;
using namespace mrcpp;

namespace mw_filter {

TEST_CASE("Interpolating filters", "[mw_filter]") {
    int maxOrder = 40;
    getInterpolatingFilterCache(ifilters);

    SECTION("Unit matrix") {
        double off_diag_row = 0.0;
        double off_diag_col = 0.0;
        for (int k = 1; k < maxOrder; k++) {
            for (int i = 0; i < k+1; i++) {
                for (int j = i; j < k+1; j++) {
                    const MatrixXd &F = ifilters.get(k).getFilter();
                    double sc = F.col(i).dot(F.col(j));
                    double sr = F.row(i).dot(F.row(j));
                    if (i == j) {
                        REQUIRE( sc == Approx(1.0) );
                        REQUIRE( sr == Approx(1.0) );
                    } else {
                        off_diag_col += std::abs(sc);
                        off_diag_row += std::abs(sr);
                    }
                }
            }
        }
        REQUIRE( off_diag_col == Approx(0.0).margin(1.0e-12) );
        REQUIRE( off_diag_row == Approx(0.0).margin(1.0e-12) );
    }
}

} // namespace
