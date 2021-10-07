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

/**
 *
 *  Base class for general polynomials with reasonably advanced
 * properties. The Polynomial class(es) are not implemented in the
 * most efficient manner, because they are only evaluated a fixed
 * number of times in a few predefined points, and all other
 * evaluations are done by linear transformations. PolynomialCache
 * implements the fast, and static const versions of the various
 * 4Polynomials.
 */

#pragma once

#include <vector>

#include <Eigen/Core>

#include "RepresentableFunction.h"

namespace mrcpp {

class Polynomial : public RepresentableFunction<1> {
public:
    Polynomial(int k = 0, const double *a = nullptr, const double *b = nullptr);
    Polynomial(int k, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(k, a.data(), b.data()) {}
    Polynomial(const Eigen::VectorXd &c, const double *a = nullptr, const double *b = nullptr);
    Polynomial(const Eigen::VectorXd &c, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(c, a.data(), b.data()) {}
    Polynomial(double c, int k = 0, const double *a = nullptr, const double *b = nullptr);
    Polynomial(double c, int k, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(c, k, a.data(), b.data()) {}
    Polynomial(const Polynomial &poly);
    Polynomial &operator=(const Polynomial &poly);
    virtual ~Polynomial() = default;

    double evalf(double x) const;
    double evalf(const Coord<1> &r) const { return evalf(r[0]); }

    double getScaledLowerBound() const;
    double getScaledUpperBound() const;

    void normalize();
    double calcSquareNorm();

    double getTranslation() const { return this->L; }
    double getDilation() const { return this->N; }

    void setDilation(double n) { this->N = n; }
    void setTranslation(double l) { this->L = l; }
    void dilate(double n) { this->N *= n; }
    void translate(double l) { this->L += this->N*l; }

    int size() const { return this->coefs.size(); } ///< Length of coefs vector
    int getOrder() const;
    void clearCoefs() { this->coefs = Eigen::VectorXd::Zero(1); }
    void setZero() { this->coefs = Eigen::VectorXd::Zero(this->coefs.size()); }
    void setCoefs(const Eigen::VectorXd &c) { this->coefs = c; }

    Eigen::VectorXd &getCoefs() { return this->coefs; }
    const Eigen::VectorXd &getCoefs() const { return this->coefs; }

    Polynomial calcDerivative() const;
    Polynomial calcAntiDerivative() const;

    void calcDerivativeInPlace();
    void calcAntiDerivativeInPlace();

    double integrate(const double *a = 0, const double *b = 0) const;
    double innerProduct(const Polynomial &p) const;

    void addInPlace(double c, const Polynomial &Q);
    Polynomial add(double c, const Polynomial &Q) const;

    Polynomial operator*(double c) const;
    Polynomial operator*(const Polynomial &Q) const;
    Polynomial operator+(const Polynomial &Q) const { return add(1.0, Q); }
    Polynomial operator-(const Polynomial &Q) const { return add(-1.0, Q); }
    Polynomial &operator*=(double c);
    Polynomial &operator*=(const Polynomial &Q);
    Polynomial &operator+=(const Polynomial &Q);
    Polynomial &operator-=(const Polynomial &Q);

protected:
    double N;              ///< Dilation coeff
    double L;              ///< Translation coeff
    Eigen::VectorXd coefs; ///< Expansion coefficients
};

} // namespace mrcpp
