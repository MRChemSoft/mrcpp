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

#include "functions/Polynomial.h"

using namespace Eigen;
using namespace mrcpp;

namespace polynomial {

TEST_CASE("Polynomial constructors", "[poly_constructor], [polynomials]") {
    double a = 0.0;
    double b = 2.0;
    Vector4d c = {0.0, 1.0, 0.0, 0.0};
    Polynomial P(c, &a, &b);

    SECTION("Constructor") {
        REQUIRE(P.getOrder() == 1);
        REQUIRE(P.getCoefs()[0] == Catch::Approx(0.0));
        REQUIRE(P.getDilation() == Catch::Approx(1.0));
        REQUIRE(P.getTranslation() == Catch::Approx(0.0));
        REQUIRE(P.isBounded());
    }

    SECTION("Copy constructor") {
        Polynomial Q(P);
        REQUIRE(Q.getOrder() == 1);
        REQUIRE(Q.getCoefs()[0] == Catch::Approx(0.0));
        REQUIRE(Q.getCoefs()[1] == Catch::Approx(1.0));
        REQUIRE(Q.getLowerBound(0) == Catch::Approx(0.0));
        REQUIRE(Q.getUpperBound(0) == Catch::Approx(2.0));
    }

    SECTION("Default constructor") {
        Polynomial Q;
        REQUIRE(Q.getOrder() == 0);
        REQUIRE(Q.getCoefs()[0] == Catch::Approx(0.0));
        REQUIRE(Q.getDilation() == Catch::Approx(1.0));
        REQUIRE(Q.getTranslation() == Catch::Approx(0.0));
        REQUIRE_FALSE(Q.isBounded());

        SECTION("Assignment operator") {
            Q = P;
            REQUIRE(Q.getOrder() == 1);
            REQUIRE(Q.getCoefs()[0] == Catch::Approx(0.0));
            REQUIRE(Q.getCoefs()[1] == Catch::Approx(1.0));
            REQUIRE_FALSE(Q.isBounded());
        }
    }
}

TEST_CASE("Polynomial evaluation", "[poly_evalf], [polynomials]") {
    double a = 0.0;
    double b = 2.0;
    Vector3d c = {1.0, 0.0, 2.0};
    Polynomial P(c, &a, &b);

    SECTION("Evaluation within bounds") {
        Coord<1> x{1.5};
        double calc_val = P.evalf(x);
        double ref_val = 1.0 + 2.0 * x[0] * x[0];
        REQUIRE(calc_val == Catch::Approx(ref_val));
    }

    SECTION("Evaluation out of bounds") {
        Coord<1> x{2.5};
        double calc_val = P.evalf(x);
        double ref_val = 0.0;
        REQUIRE(calc_val == Catch::Approx(ref_val));
    }
}

SCENARIO("Polynomials can be scaled and translated", "[poly_scale], [polynomials]") {
    GIVEN("A bounded polynomial P") {
        double a = -1.0;
        double b = 1.0;
        Vector3d c = {0.0, 1.0, 1.0};
        Polynomial P(c, &a, &b);
        WHEN("P is translated then dilated") {
            double n = 2.0;
            double l = 1.0;
            P.translate(l);
            P.dilate(n);
            THEN("The dilation changes") { REQUIRE(P.getDilation() == Catch::Approx(2.0)); }
            THEN("The translation changes") { REQUIRE(P.getTranslation() == Catch::Approx(1.0)); }
            THEN("The scaled bounds change") {
                REQUIRE(P.getScaledLowerBound() == Catch::Approx(0.0));
                REQUIRE(P.getScaledUpperBound() == Catch::Approx(1.0));
            }
            THEN("The unscaled bounds don't change") {
                REQUIRE(P.getLowerBound(0) == Catch::Approx(-1.0));
                REQUIRE(P.getUpperBound(0) == Catch::Approx(1.0));
            }
            THEN("The scaled evaluation is known") {
                Coord<1> x{0.3};
                double calc_val = P.evalf(x);
                double ref_val = (2.0 * x[0] - 1.0) + (2.0 * x[0] - 1.0) * (2.0 * x[0] - 1.0);
                REQUIRE(calc_val == Catch::Approx(ref_val));
            }
        }
        WHEN("P is dilated then translated") {
            double n = 2.0;
            double l = 1.0;
            P.dilate(n);
            P.translate(l);
            THEN("The dilation changes") { REQUIRE(P.getDilation() == Catch::Approx(2.0)); }
            THEN("The translation changes") { REQUIRE(P.getTranslation() == Catch::Approx(2.0)); }
            THEN("The scaled bounds change") {
                REQUIRE(P.getScaledLowerBound() == Catch::Approx(0.5));
                REQUIRE(P.getScaledUpperBound() == Catch::Approx(1.5));
            }
            THEN("The unscaled bounds don't change") {
                REQUIRE(P.getLowerBound(0) == Catch::Approx(-1.0));
                REQUIRE(P.getUpperBound(0) == Catch::Approx(1.0));
            }
            THEN("The scaled evaluation is known") {
                Coord<1> x{0.7};
                double calc_val = P.evalf(x);
                double ref_val = (2.0 * x[0] - 2.0) + (2.0 * x[0] - 2.0) * (2.0 * x[0] - 2.0);
                REQUIRE(calc_val == Catch::Approx(ref_val));
            }
        }
    }
}

SCENARIO("Polynomials can be added and multiplied", "[poly_arithmetics], [polynomials]") {
    GIVEN("Two polynomials P and Q") {
        Vector3d pCoefs = {0.0, 1.0, 1.0};
        Vector2d qCoefs = {0.0, 1.0};
        Polynomial P(pCoefs);
        Polynomial Q(qCoefs);

        WHEN("Q += P") {
            Q += P;
            THEN("The coefficients of P are unchanged") {
                REQUIRE(P.getOrder() == 2);
                REQUIRE(P.getCoefs()[0] == Catch::Approx(0.0));
                REQUIRE(P.getCoefs()[1] == Catch::Approx(1.0));
                REQUIRE(P.getCoefs()[2] == Catch::Approx(1.0));
            }
            THEN("The coefficients of Q change") {
                REQUIRE(Q.getOrder() == 2);
                REQUIRE(Q.getCoefs()[0] == Catch::Approx(0.0));
                REQUIRE(Q.getCoefs()[1] == Catch::Approx(2.0));
                REQUIRE(Q.getCoefs()[2] == Catch::Approx(1.0));
            }
        }

        WHEN("R = P + Q") {
            Polynomial R;
            R = P + Q;
            THEN("The coefficients of R are known") {
                REQUIRE(R.getOrder() == 2);
                REQUIRE(R.getCoefs()[0] == Catch::Approx(0.0));
                REQUIRE(R.getCoefs()[1] == Catch::Approx(2.0));
                REQUIRE(R.getCoefs()[2] == Catch::Approx(1.0));
            }
        }

        WHEN("P *= 2.0") {
            P *= 2.0;
            THEN("The coefficients of P change") {
                REQUIRE(P.getOrder() == 2);
                REQUIRE(P.getCoefs()[0] == Catch::Approx(0.0));
                REQUIRE(P.getCoefs()[1] == Catch::Approx(2.0));
                REQUIRE(P.getCoefs()[2] == Catch::Approx(2.0));
            }
        }

        WHEN("R = P*3.0") {
            Polynomial R;
            R = P * 3.0;
            THEN("The coefficients of R are known") {
                REQUIRE(R.getOrder() == 2);
                REQUIRE(R.getCoefs()[0] == Catch::Approx(0.0));
                REQUIRE(R.getCoefs()[1] == Catch::Approx(3.0));
                REQUIRE(R.getCoefs()[2] == Catch::Approx(3.0));
            }
        }

        WHEN("Q *= P") {
            Q *= P;
            THEN("The coefficients of P are unchanged") {
                REQUIRE(P.getOrder() == 2);
                REQUIRE(P.getCoefs()[0] == Catch::Approx(0.0));
                REQUIRE(P.getCoefs()[1] == Catch::Approx(1.0));
                REQUIRE(P.getCoefs()[2] == Catch::Approx(1.0));
            }
            THEN("The coefficients of Q change") {
                REQUIRE(Q.getOrder() == 3);
                REQUIRE(Q.getCoefs()[0] == Catch::Approx(0.0));
                REQUIRE(Q.getCoefs()[1] == Catch::Approx(0.0));
                REQUIRE(Q.getCoefs()[2] == Catch::Approx(1.0));
                REQUIRE(Q.getCoefs()[3] == Catch::Approx(1.0));
            }
        }

        WHEN("R = P * Q") {
            Polynomial R;
            R = P * Q;
            THEN("The coefficients of R are known") {
                REQUIRE(R.getCoefs()[0] == Catch::Approx(0.0));
                REQUIRE(R.getCoefs()[1] == Catch::Approx(0.0));
                REQUIRE(R.getCoefs()[2] == Catch::Approx(1.0));
                REQUIRE(R.getCoefs()[3] == Catch::Approx(1.0));
            }
        }
    }
}

TEST_CASE("Polynomial differentiation", "[poly_diff], [polynomials]") {
    Vector3d c = {0.0, 1.0, 2.0};
    Polynomial P(c);

    SECTION("Derivative in place") {
        P.calcDerivativeInPlace();
        REQUIRE(P.getOrder() == 1);
        REQUIRE(P.getCoefs()[0] == Catch::Approx(1.0));
        REQUIRE(P.getCoefs()[1] == Catch::Approx(4.0));
    }

    SECTION("Derivative") {
        Polynomial Q = P.calcDerivative();
        REQUIRE(Q.getOrder() == 1);
        REQUIRE(Q.getCoefs()[0] == Catch::Approx(1.0));
        REQUIRE(Q.getCoefs()[1] == Catch::Approx(4.0));
    }
}

TEST_CASE("Polynomial integration", "[poly_int], [polynomials]") {
    Vector3d c = {0.0, 1.0, 2.0};
    Polynomial P(c);

    SECTION("Antiderivative in place") {
        P.calcAntiDerivativeInPlace();
        REQUIRE(P.getOrder() == 3);
        REQUIRE(P.getCoefs()[0] == Catch::Approx(0.0));
        REQUIRE(P.getCoefs()[1] == Catch::Approx(0.0));
        REQUIRE(P.getCoefs()[2] == Catch::Approx(1.0 / 2.0));
        REQUIRE(P.getCoefs()[3] == Catch::Approx(2.0 / 3.0));
    }

    SECTION("Antiderivative") {
        Polynomial Q = P.calcAntiDerivative();
        REQUIRE(Q.getOrder() == 3);
        REQUIRE(Q.getCoefs()[0] == Catch::Approx(0.0));
        REQUIRE(Q.getCoefs()[1] == Catch::Approx(0.0));
        REQUIRE(Q.getCoefs()[2] == Catch::Approx(1.0 / 2.0));
        REQUIRE(Q.getCoefs()[3] == Catch::Approx(2.0 / 3.0));
    }

    GIVEN("A bounded polynomial P on [-1.0, 1.0]") {
        double a = -1.0;
        double b = 1.0;
        P.setBounds(&a, &b);
        THEN("P can be integrated on its full domain") {
            double calc_int = P.integrate();
            double ref_int = 4.0 / 3.0;
            REQUIRE(calc_int == Catch::Approx(ref_int));
        }
        THEN("P can be integrated on the subdomain [0.0, 1.0]") {
            a = 0.0;
            double calc_int = P.integrate(&a, &b);
            double ref_int = 7.0 / 6.0;
            REQUIRE(calc_int == Catch::Approx(ref_int));
        }
    }
}

SCENARIO("Bounded polynomials have inner products and norms", "[poly_norm], [polynomials]") {
    GIVEN("Two unbounded polynomials P and Q") {
        Vector3d c1 = {0.0, 1.0, 1.0};
        Vector2d c2 = {0.0, 1.0};
        Polynomial P(c1);
        Polynomial Q(c2);
        THEN("The norm of P is undefined") { REQUIRE(P.calcSquareNorm() < 0.0); }
        WHEN("P is bounded") {
            double a = -1.0;
            double b = 1.0;
            P.setBounds(&a, &b);
            THEN("The inner product <P|Q> is defined") {
                double ref_inner = 2.0 / 3.0;
                double calc_inner = P.innerProduct(Q);
                REQUIRE(calc_inner == Catch::Approx(ref_inner));
            }
            THEN("The norm of P is defined") {
                double ref_norm = 16.0 / 15.0;
                double calc_norm = P.calcSquareNorm();
                REQUIRE(calc_norm == Catch::Approx(ref_norm));
            }
            THEN("P can be normalized") {
                P.normalize();
                double calc_norm = P.calcSquareNorm();
                REQUIRE(calc_norm == Catch::Approx(1.0));
            }
        }
    }
}

} // namespace polynomial
