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

#include "functions/GaussExp.h"
#include "functions/GaussFunc.h"

#include "utils/details.h"

using namespace mrcpp;

namespace periodify_gaussians {

template <int D> void testIsNarrowGaussiansMadePeriodic();
template <int D> void testIsWideGaussiansMadePeriodic();
template <int D> void testIsGaussExpMadePeriodic();

template <int D> bool isPeriodic(const Gaussian<D> &func) {
    if (func.getPeriod() == std::array<double, D>{})
        return false;
    else
        return true;
}

SCENARIO("Given a narrow Gaussian is it periodic?", "[periodic_narrow_gaussian]") {
    GIVEN("Narrow Gaussian in 1D") { testIsNarrowGaussiansMadePeriodic<1>(); }
    GIVEN("Narrow Gaussian in 2D") { testIsNarrowGaussiansMadePeriodic<2>(); }
    GIVEN("Narrow Gaussian in 3D") { testIsNarrowGaussiansMadePeriodic<3>(); }
}

SCENARIO("Given a wide Gaussian is it periodic?", "[periodic_wide_gaussian]") {
    GIVEN("Wide Gaussian in 1D") { testIsWideGaussiansMadePeriodic<1>(); }
    GIVEN("Wide Gaussian in 2D") { testIsWideGaussiansMadePeriodic<2>(); }
    GIVEN("Wide Gaussian in 3D") { testIsWideGaussiansMadePeriodic<3>(); }
}

SCENARIO("Given a GaussExp, can the Gaussians be made periodic?", "[periodic_gaussExp]") {
    GIVEN("GaussExp in 1D") { testIsGaussExpMadePeriodic<1>(); }
    GIVEN("GaussExp in 2D") { testIsGaussExpMadePeriodic<2>(); }
    GIVEN("GaussExp in 3D") { testIsGaussExpMadePeriodic<3>(); }
}

template <int D> void testIsNarrowGaussiansMadePeriodic() {

    auto beta = 1.0e4;
    auto alpha = std::pow(beta / mrcpp::pi, D / 2.0);
    double pos_data[3] = {-0.2, 0.5, 1.0};

    auto period = std::array<double, D>{};
    period.fill(2.0);

    auto pos = mrcpp::details::convert_to_std_array<double, D>(pos_data);
    auto new_pos = pos;
    new_pos[0] += period[0]; // Zero is selected such that it will work inn 1, 2 and 3D

    WHEN("The gaussian is generated it is constructed, it is not periodic") {
        auto gauss = GaussFunc<D>(beta, alpha, pos);
        REQUIRE(not isPeriodic<D>(gauss));
        THEN("The Gaussian should be close to zero at the bondary of the unit cell") {
            REQUIRE(gauss.evalf(new_pos) == Approx(0.0));
        }
        AND_WHEN("The gaussian is made periodc, it is periodic") {
            gauss.makePeriodic(period);
            REQUIRE(isPeriodic<D>(gauss));
            THEN("The gaussian should have the same value at -0.2 and 1.8") {
                REQUIRE(gauss.evalf(pos) == gauss.evalf(new_pos));
            }
        }
    }
}

template <int D> void testIsWideGaussiansMadePeriodic() {

    auto beta = 1.0;
    auto alpha = std::pow(beta / mrcpp::pi, D / 2.0);
    double pos_data[3] = {-0.2, 0.5, 1.0};

    auto period = std::array<double, D>{};
    period.fill(2.0);

    auto pos = mrcpp::details::convert_to_std_array<double, D>(pos_data);
    auto new_pos = pos;
    new_pos[0] += period[0]; // Zero is selected such that it will work inn 1, 2 and 3D

    WHEN("The gaussian is generated it is constructed, it is not periodic") {
        auto gauss = GaussFunc<D>(beta, alpha, pos);
        REQUIRE(not isPeriodic<D>(gauss));
        THEN("The Gaussian should be close to zero at the bondary of the unit cell") {
            REQUIRE(gauss.evalf(new_pos) != Approx(0.0));
        }
        AND_WHEN("The gaussian is made periodc, it is periodic") {
            gauss.makePeriodic(period);
            REQUIRE(isPeriodic<D>(gauss));
            THEN("The gaussian should have the same value at -0.2 and 1.8") {
                REQUIRE(gauss.evalf(pos) == Approx(gauss.evalf(new_pos)));
            }
        }
    }
}

template <int D> void testIsGaussExpMadePeriodic() {

    auto beta = 1.0e4;
    auto alpha = std::pow(beta / mrcpp::pi, D / 2.0);
    double pos_1_data[3] = {-0.2, 0.5, 1.0};
    double pos_2_data[3] = {-0.1, 0.3, 1.2};

    auto period = std::array<double, D>{};
    period.fill(2.0);

    auto pos_1 = mrcpp::details::convert_to_std_array<double, D>(pos_1_data);
    auto pos_2 = mrcpp::details::convert_to_std_array<double, D>(pos_2_data);
    auto gauss_1 = GaussFunc<D>(beta, alpha, pos_1);
    auto gauss_2 = GaussFunc<D>(beta, alpha, pos_2);

    auto gauss_exp = GaussExp<D>();
    WHEN("The gaussians are added they are non-periodic") {
        gauss_exp.append(gauss_1);
        gauss_exp.append(gauss_2);
        for (auto &func : gauss_exp) REQUIRE(not isPeriodic<D>(*func));

        AND_WHEN("The gaussians are make periodic") {
            gauss_exp.makePeriodic(period);
            for (auto &func : gauss_exp) REQUIRE(isPeriodic<D>(*func));
        }
    }
}

} // namespace periodify_gaussians
