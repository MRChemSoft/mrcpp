#pragma once

#include "Polynomial.h"

namespace mrcpp {

class LegendrePoly final : public Polynomial {
public:
    LegendrePoly(int k, double n = 1.0, double l = 0.0);

    Eigen::Vector2d firstDerivative(double x) const;
    Eigen::Vector3d secondDerivative(double x) const;

private:
    void computeLegendrePolynomial(int k);
};

}
