#include "catch.hpp"

#include "Polynomial.h"

using namespace Eigen;

namespace polynomial {

TEST_CASE("Polynomial constructors", "[poly_constructor], [polynomials]") {
    double a = 0.0;
    double b = 2.0;
    Vector4d c = {0.0, 1.0, 0.0, 0.0};
    Polynomial P(c, &a, &b);

    SECTION("Constructor") {
        REQUIRE( (P.getOrder() == 1) );
        REQUIRE( (P.getCoefs()[0] == Approx(0.0)) );
        REQUIRE( (P.getDilation() == Approx(1.0)) );
        REQUIRE( (P.getTranslation() == Approx(0.0)) );
        REQUIRE( P.isBounded() );
    }

    SECTION("Copy constructor") {
        Polynomial Q(P);
        REQUIRE( (Q.getOrder() == 1) );
        REQUIRE( (Q.getCoefs()[0] == Approx(0.0)) );
        REQUIRE( (Q.getCoefs()[1] == Approx(1.0)) );
        REQUIRE( (Q.getLowerBound(0) == 0.0) );
        REQUIRE( (Q.getUpperBound(0) == 2.0) );
    }

    SECTION("Default constructor") {
        Polynomial Q;
        REQUIRE( (Q.getOrder() == 0) );
        REQUIRE( (Q.getCoefs()[0] == Approx(0.0)) );
        REQUIRE( (Q.getDilation() == Approx(1.0)) );
        REQUIRE( (Q.getTranslation() == Approx(0.0)) );
        REQUIRE_FALSE( Q.isBounded() );

        SECTION("Assignment operator") {
            Q = P;
            REQUIRE( (Q.getOrder() == 1) );
            REQUIRE( (Q.getCoefs()[0] == Approx(0.0)) );
            REQUIRE( (Q.getCoefs()[1] == Approx(1.0)) );
            REQUIRE_FALSE( Q.isBounded() );
        }
    }
}

TEST_CASE("Polynomial evaluation", "[poly_evalf], [polynomials]") {
    double a = 0.0;
    double b = 2.0;
    Vector3d c = {1.0, 0.0, 2.0};
    Polynomial P(c, &a, &b);

    SECTION("Evaluation within bounds") {
        double x = 1.5;
        double calc_val = P.evalf(&x);
        double ref_val = 1.0 + 2.0*x*x;
        REQUIRE( (calc_val == Approx(ref_val)) );
    }

    SECTION("Evaluation out of bounds") {
        double x = 2.5;
        double calc_val = P.evalf(&x);
        double ref_val = 0.0;
        REQUIRE( (calc_val == Approx(ref_val)) );
    }
}

SCENARIO("Polynomials can be scaled and translated", "[poly_scale], [polynomials]") {
    GIVEN("A bounded polynomial P") {
        double a = -1.0;
        double b = 1.0;
        Vector3d c = {0.0, 1.0, 1.0};
        Polynomial P(c, &a, &b);
        WHEN("P is rescaled") {
            double n = 2.0;
            double l = 1.0;
            P.rescale(n, l);
            THEN("The dilation changes") {
                REQUIRE( (P.getDilation() == Approx(2.0)) );
            }
            THEN("The translation changes") {
                REQUIRE( (P.getTranslation() == Approx(1.0)) );
            }
            THEN("The scaled bounds change") {
                REQUIRE( (P.getScaledLowerBound() == Approx(0.0)) );
                REQUIRE( (P.getScaledUpperBound() == Approx(1.0)) );
            }
            THEN("The unscaled bounds don't change") {
                REQUIRE( (P.getLowerBound(0) == Approx(-1.0)) );
                REQUIRE( (P.getUpperBound(0) == Approx(1.0)) );
            }
            THEN("The scaled evaluation is known") {
                double x = 0.3;
                double calc_val = P.evalf(&x);
                double ref_val = (2.0*x - 1.0) + (2.0*x - 1.0)*(2.0*x - 1.0);
                REQUIRE( (calc_val == Approx(ref_val)) );
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
                REQUIRE( (P.getOrder() == 2) );
                REQUIRE( (P.getCoefs()[0] == Approx(0.0)) );
                REQUIRE( (P.getCoefs()[1] == Approx(1.0)) );
                REQUIRE( (P.getCoefs()[2] == Approx(1.0)) );
            }
            THEN("The coefficients of Q change") {
                REQUIRE( (Q.getOrder() == 2) );
                REQUIRE( (Q.getCoefs()[0] == Approx(0.0)) );
                REQUIRE( (Q.getCoefs()[1] == Approx(2.0)) );
                REQUIRE( (Q.getCoefs()[2] == Approx(1.0)) );
            }
        }

        WHEN("R = P + Q") {
            Polynomial R;
            R = P+Q;
            THEN("The coefficients of R are known") {
                REQUIRE( (R.getOrder() == 2) );
                REQUIRE( (R.getCoefs()[0] == Approx(0.0)) );
                REQUIRE( (R.getCoefs()[1] == Approx(2.0)) );
                REQUIRE( (R.getCoefs()[2] == Approx(1.0)) );
            }
        }

        WHEN("P *= 2.0") {
            P *= 2.0;
            THEN("The coefficients of P change") {
                REQUIRE( (P.getOrder() == 2) );
                REQUIRE( (P.getCoefs()[0] == Approx(0.0)) );
                REQUIRE( (P.getCoefs()[1] == Approx(2.0)) );
                REQUIRE( (P.getCoefs()[2] == Approx(2.0)) );
            }
        }

        WHEN("R = P*3.0") {
            Polynomial R;
            R = P*3.0;
            THEN("The coefficients of R are known") {
                REQUIRE( (R.getOrder() == 2) );
                REQUIRE( (R.getCoefs()[0] == Approx(0.0)) );
                REQUIRE( (R.getCoefs()[1] == Approx(3.0)) );
                REQUIRE( (R.getCoefs()[2] == Approx(3.0)) );
            }
        }

        WHEN("Q *= P") {
            Q *= P;
            THEN("The coefficients of P are unchanged") {
                REQUIRE( (P.getOrder() == 2) );
                REQUIRE( (P.getCoefs()[0] == Approx(0.0)) );
                REQUIRE( (P.getCoefs()[1] == Approx(1.0)) );
                REQUIRE( (P.getCoefs()[2] == Approx(1.0)) );
            }
            THEN("The coefficients of Q change") {
                REQUIRE( (Q.getOrder() == 3) );
                REQUIRE( (Q.getCoefs()[0] == Approx(0.0)) );
                REQUIRE( (Q.getCoefs()[1] == Approx(0.0)) );
                REQUIRE( (Q.getCoefs()[2] == Approx(1.0)) );
                REQUIRE( (Q.getCoefs()[3] == Approx(1.0)) );
            }
        }

        WHEN("R = P * Q") {
            Polynomial R;
            R = P*Q;
            THEN("The coefficients of R are known") {
                VectorXd &cr = R.getCoefs();
                REQUIRE( (R.getCoefs()[0] == Approx(0.0)) );
                REQUIRE( (R.getCoefs()[1] == Approx(0.0)) );
                REQUIRE( (R.getCoefs()[2] == Approx(1.0)) );
                REQUIRE( (R.getCoefs()[3] == Approx(1.0)) );
            }
        }
    }
}

TEST_CASE("Polynomial differentiation", "[poly_diff], [polynomials]") {
    Vector3d c = {0.0, 1.0, 2.0};
    Polynomial P(c);

    SECTION("Derivative in place") {
        P.calcDerivativeInPlace();
        REQUIRE( (P.getOrder() == 1) );
        REQUIRE( (P.getCoefs()[0] == Approx(1.0)) );
        REQUIRE( (P.getCoefs()[1] == Approx(4.0)) );
    }

    SECTION("Derivative") {
        Polynomial Q = P.calcDerivative();
        REQUIRE( (Q.getOrder() == 1) );
        REQUIRE( (Q.getCoefs()[0] == Approx(1.0)) );
        REQUIRE( (Q.getCoefs()[1] == Approx(4.0)) );
    }
}

TEST_CASE("Polynomial integration", "[poly_int], [polynomials]") {
    Vector3d c = {0.0, 1.0, 2.0};
    Polynomial P(c);

    SECTION("Antiderivative in place") {
        P.calcAntiDerivativeInPlace();
        REQUIRE( (P.getOrder() == 3) );
        REQUIRE( (P.getCoefs()[0] == Approx(0.0)) );
        REQUIRE( (P.getCoefs()[1] == Approx(0.0)) );
        REQUIRE( (P.getCoefs()[2] == Approx(1.0/2.0)) );
        REQUIRE( (P.getCoefs()[3] == Approx(2.0/3.0)) );
    }

    SECTION("Antiderivative") {
        Polynomial Q = P.calcAntiDerivative();
        REQUIRE( (Q.getOrder() == 3) );
        REQUIRE( (Q.getCoefs()[0] == Approx(0.0)) );
        REQUIRE( (Q.getCoefs()[1] == Approx(0.0)) );
        REQUIRE( (Q.getCoefs()[2] == Approx(1.0/2.0)) );
        REQUIRE( (Q.getCoefs()[3] == Approx(2.0/3.0)) );
    }

    GIVEN("A bounded polynomial P on [-1.0, 1.0]") {
        double a = -1.0;
        double b = 1.0;
        P.setBounds(&a, &b);
        THEN("P can be integrated on its full domain") {
            double calc_int = P.integrate();
            double ref_int = 4.0/3.0;
            REQUIRE( (calc_int == Approx(ref_int)) );
        }
        THEN("P can be integrated on the subdomain [0.0, 1.0]") {
            a = 0.0;
            double calc_int = P.integrate(&a, &b);
            double ref_int = 7.0/6.0;
            REQUIRE( (calc_int == Approx(ref_int)) );
        }
    }
}

SCENARIO("Bounded polynomials have inner products and norms", "[poly_norm], [polynomials]") {
    GIVEN("Two unbounded polynomials P and Q") {
        Vector3d c1 = {0.0, 1.0, 1.0};
        Vector2d c2 = {0.0, 1.0};
        Polynomial P(c1);
        Polynomial Q(c2);
        THEN("The norm of P is undefined") {
            REQUIRE( (P.calcSquareNorm() < 0.0) );
        }
        WHEN("P is bounded") {
            double a = -1.0;
            double b = 1.0;
            P.setBounds(&a, &b);
            THEN("The inner product <P|Q> is defined") {
                double ref_inner = 2.0/3.0;
                double calc_inner = P.innerProduct(Q);
                REQUIRE( (calc_inner == Approx(ref_inner)) );
            }
            THEN("The norm of P is defined") {
                double ref_norm = 16.0/15.0;
                double calc_norm = P.calcSquareNorm();
                REQUIRE( (calc_norm == Approx(ref_norm)) );
            }
            THEN("P can be normalized") {
                P.normalize();
                double calc_norm = P.calcSquareNorm();
                REQUIRE( (calc_norm == Approx(1.0)) );
            }
        }
    }
}

} // namespace
