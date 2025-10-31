/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3 of the License, or
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

#include <cfloat>

#include "Polynomial.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

Polynomial::Polynomial(int k, const double *a, const double *b)
        : RepresentableFunction<1, double>(a, b) {
    assert(k >= 0);
    this->N = 1.0;
    this->L = 0.0;
    this->coefs = VectorXd::Zero(k + 1);
}

Polynomial::Polynomial(double c, int k, const double *a, const double *b)
        : RepresentableFunction<1>(a, b) {
    this->N = 1.0;
    this->L = 0.0;
    this->coefs = math_utils::get_binomial_coefs(k);
    for (int i = 0; i <= k; i++) { this->coefs[i] *= std::pow(c, k - i); }
}

Polynomial::Polynomial(const VectorXd &c, const double *a, const double *b)
        : RepresentableFunction<1>(a, b) {
    this->N = 1.0;
    this->L = 0.0;
    setCoefs(c);
}

Polynomial::Polynomial(const Polynomial &poly)
        : RepresentableFunction<1>(poly) {
    this->N = poly.N;
    this->L = poly.L;
    this->coefs = poly.coefs;
}

Polynomial &Polynomial::operator=(const Polynomial &poly) {
    RepresentableFunction<1>::operator=(poly);
    this->N = poly.N;
    this->L = poly.L;
    this->coefs = poly.coefs;
    return *this;
}

double Polynomial::evalf(double x) const {
    if (isBounded()) {
        if (x < this->getScaledLowerBound()) return 0.0;
        if (x > this->getScaledUpperBound()) return 0.0;
    }
    double xp = 1.0;
    double y = 0.0;
    for (int k = 0; k < getOrder() + 1; k++) {
        y += (xp * this->coefs[k]);
        xp *= this->N * x - this->L;
    }
    return y;
}

double Polynomial::getScaledLowerBound() const {
    if (not isBounded()) MSG_ERROR("Unbounded polynomial");
    return (1.0 / this->N * (this->A[0] + this->L));
}

double Polynomial::getScaledUpperBound() const {
    if (not isBounded()) MSG_ERROR("Unbounded polynomial");
    return (1.0 / this->N * (this->B[0] + this->L));
}

void Polynomial::normalize() {
    double sqNorm = calcSquareNorm();
    if (sqNorm < 0.0) MSG_ABORT("Cannot normalize polynomial");
    (*this) *= 1.0 / std::sqrt(sqNorm);
}

double Polynomial::calcSquareNorm() {
    double sqNorm = -1.0;
    if (isBounded()) { sqNorm = this->innerProduct(*this); }
    return sqNorm;
}

int Polynomial::getOrder() const {
    int n = 0;
    for (int i = 0; i < this->coefs.size(); i++) {
        if (std::abs(this->coefs[i]) > MachineZero) { n = i; }
    }
    return n;
}

Polynomial &Polynomial::operator*=(double c) {
    this->coefs = c * this->coefs;
    return *this;
}

Polynomial &Polynomial::operator*=(const Polynomial &Q) {
    Polynomial &P = *this;
    if (std::abs(P.getDilation() - Q.getDilation()) > MachineZero) { MSG_ERROR("Polynomials not defined on same scale."); }
    if (std::abs(P.getTranslation() - Q.getTranslation()) > MachineZero) { MSG_ERROR("Polynomials not defined on same translation."); }

    int P_order = P.getOrder();
    int Q_order = Q.getOrder();
    int new_order = P_order + Q_order;
    VectorXd coefs = VectorXd::Zero(new_order + 1);
    for (int i = 0; i < P_order + 1; i++) {
        for (int j = 0; j < Q_order + 1; j++) { coefs(i + j) += P.coefs(i) * Q.coefs(j); }
    }
    P.setCoefs(coefs);
    return P;
}

Polynomial Polynomial::operator*(double c) const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q *= c;
    return Q;
}

Polynomial Polynomial::operator*(const Polynomial &Q) const {
    const Polynomial &P = *this;
    Polynomial R;
    R = P;
    R *= Q;
    return R;
}

Polynomial &Polynomial::operator+=(const Polynomial &Q) {
    this->addInPlace(1.0, Q);
    return *this;
}

Polynomial &Polynomial::operator-=(const Polynomial &Q) {
    this->addInPlace(-1.0, Q);
    return *this;
}

void Polynomial::addInPlace(double c, const Polynomial &Q) {
    Polynomial &P = *this;
    if (std::abs(P.getDilation() - Q.getDilation()) > MachineZero) { MSG_ERROR("Polynomials not defined on same scale."); }
    if (std::abs(P.getTranslation() - Q.getTranslation()) > MachineZero) { MSG_ERROR("Polynomials not defined on same translation."); }

    int P_order = P.getOrder();
    int Q_order = Q.getOrder();
    int new_order = std::max(P_order, Q_order);
    VectorXd newCoefs = VectorXd::Zero(new_order + 1);

    for (int i = 0; i < new_order + 1; i++) {
        if (i <= P_order) { newCoefs[i] += P.getCoefs()[i]; }
        if (i <= Q_order) { newCoefs[i] += c * Q.getCoefs()[i]; }
    }
    P.setCoefs(newCoefs);
}

Polynomial Polynomial::add(double c, const Polynomial &Q) const {
    const Polynomial &P = *this;
    Polynomial R;
    R = P;
    R.addInPlace(c, Q);
    return R;
}

Polynomial Polynomial::calcDerivative() const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q.calcDerivativeInPlace();
    return Q;
}

void Polynomial::calcDerivativeInPlace() {
    Polynomial &P = *this;
    int P_order = P.getOrder();
    const VectorXd &oldCoefs = P.getCoefs();
    VectorXd newCoefs = VectorXd::Zero(P_order);
    for (int i = 0; i < newCoefs.size(); i++) { newCoefs[i] = double(i + 1) * oldCoefs[i + 1]; }
    P.setCoefs(newCoefs);
}

Polynomial Polynomial::calcAntiDerivative() const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q.calcAntiDerivativeInPlace();
    return Q;
}

void Polynomial::calcAntiDerivativeInPlace() {
    Polynomial &P = *this;
    int P_order = P.getOrder();
    const VectorXd &oldCoefs = P.getCoefs();
    VectorXd newCoefs = VectorXd::Zero(P_order + 2);
    newCoefs[0] = 0.0;
    newCoefs[1] = oldCoefs[0];
    for (int i = 2; i < newCoefs.size(); i++) { newCoefs[i] = 1.0 / i * oldCoefs[i - 1]; }
    P.setCoefs(newCoefs);
}

double Polynomial::integrate(const double *a, const double *b) const {
    double lb = -DBL_MAX, ub = DBL_MAX;
    if (this->isBounded()) {
        lb = getScaledLowerBound();
        ub = getScaledUpperBound();
    } else {
        if (a == nullptr) MSG_ERROR("Polynomial without bounds");
        if (b == nullptr) MSG_ERROR("Polynomial without bounds");
    }
    if (a != nullptr) lb = std::max(a[0], lb);
    if (b != nullptr) ub = std::min(b[0], ub);
    if (lb >= ub) return 0.0;

    double sfac = 1.0 / this->N;
    Polynomial antidiff = calcAntiDerivative();
    return sfac * (antidiff.evalf(ub) - antidiff.evalf(lb));
}

double Polynomial::innerProduct(const Polynomial &Q) const {
    const Polynomial &P = *this;
    if (not P.isBounded()) MSG_ERROR("Unbounded polynomial");
    Polynomial pq = P * Q;
    pq.setBounds(P.getLowerBounds(), P.getUpperBounds());
    return pq.integrate();
}

} // namespace mrcpp