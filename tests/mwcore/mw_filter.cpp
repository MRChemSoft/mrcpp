#include "catch.hpp"

#include "mwcore/FilterCache.h"
#include "constants.h"

using namespace Eigen;
using namespace mrcpp;

namespace mw_filter {

TEST_CASE("Interpolating filters", "[mw_filter]") {
    const double thrs = 100*MachineZero;

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
                        REQUIRE( std::abs(sc - 1.0) < thrs);
                        REQUIRE( std::abs(sr - 1.0) < thrs);
                    } else {
                        off_diag_col += fabs(sc);
                        off_diag_row += fabs(sr);
                    }
                }
            }
        }
        REQUIRE( std::abs(off_diag_col) < thrs );
        REQUIRE( std::abs(off_diag_row) < thrs );
    }
}

} // namespace
