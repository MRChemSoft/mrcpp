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

#include <cmath>

#include "BoysFunction.h"
#include "function_utils.h"
#include "GaussExp.h"
#include "GaussFunc.h"
#include "GaussPoly.h"
#include "Polynomial.h"
#include "utils/Printer.h"
#include "utils/details.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

template <int D> Gaussian<D> *GaussFunc<D>::copy() const {
    auto *gauss = new GaussFunc<D>(*this);
    return gauss;
}

template <int D> double GaussFunc<D>::evalf(const Coord<D> &r) const {
    if (this->getScreen()) {
        for (int d = 0; d < D; d++) {
            if (r[d] < this->A[d] or r[d] > this->B[d]) { return 0.0; }
        }
    }
    double q2 = 0.0, p2 = 1.0;
    for (int d = 0; d < D; d++) {
        double q = r[d] - this->pos[d];
        q2 += this->alpha[d] * q * q;
        if (this->power[d] == 0) {
            continue;
        } else if (this->power[d] == 1) {
            p2 *= q;
        } else {
            p2 *= std::pow(q, this->power[d]);
        }
    }
    return this->coef * p2 * std::exp(-q2);
}

template <int D> double GaussFunc<D>::evalf1D(double r, int d) const {
    if (this->getScreen()) {
        if ((r < this->A[d]) or (r > this->B[d])) { return 0.0; }
    }
    double q = (r - this->pos[d]);
    double q2 = q * q;

    double p2;
    if (this->power[d] == 0) {
        p2 = 1.0;
    } else if (this->power[d] == 1) {
        p2 = q;
    } else {
        p2 = std::pow(q, this->power[d]);
    }
    double result = p2 * std::exp(-this->alpha[d] * q2);
    if (d == 0) result *= this->coef;
    return result;
}

template <int D> double GaussFunc<D>::calcSquareNorm() const {
    double norm = 1.0;
    for (int d = 0; d < D; d++) {
        double a = 2.0 * this->alpha[d];
        double sq_norm = 1.0;
        int p = this->power[d];
        if (p > 0) {
            int i = 2 * p - 1;
            while (i > 0) {
                sq_norm = i * sq_norm / (2.0 * a);
                i = i - 2;
            }
        }
        a = pi / a;
        sq_norm *= std::sqrt(a);
        norm *= sq_norm;
    }
    return norm * this->coef * this->coef;
}

template<int D> GaussExp<D> GaussFunc<D>::asGaussExp() const {
    GaussExp<D> gexp;
    gexp.append(*this);
    return gexp;
}

template <int D> GaussPoly<D> GaussFunc<D>::differentiate(int dir) const {
    GaussPoly<D> result(*this);
    int oldPow = this->getPow(dir);

    Polynomial newPoly(oldPow + 1);
    newPoly.getCoefs()[oldPow + 1] = -2.0 * this->getExp()[dir];
    if (oldPow > 0) { newPoly.getCoefs()[oldPow - 1] = oldPow; }
    result.setPoly(dir, newPoly);
    ;
    return result;
}

template <int D> void GaussFunc<D>::multInPlace(const GaussFunc<D> &rhs) {
    GaussFunc<D> &lhs = *this;
    for (int d = 0; d < D; d++) {
        if (lhs.getPos()[d] != rhs.getPos()[d]) {
            MSG_ABORT("Cannot multiply GaussFuncs of different center in-place");
        }
    }
    double newCoef = lhs.getCoef() * rhs.getCoef();
    std::array<double, D> newExp;
    auto lhsExp = lhs.getExp();
    auto rhsExp = rhs.getExp();
    for (int d = 0; d < D; d++) { newExp[d] = lhsExp[d] + rhsExp[d]; }

    std::array<int, D> newPow;
    for (int d = 0; d < D; d++) { newPow[d] = lhs.getPow(d) + rhs.getPow(d); }
    this->setCoef(newCoef);
    this->setExp(newExp);
    this->setPow(newPow);
}

/** @brief Multiply two GaussFuncs
 *  @param[in] this: Left hand side of multiply
 *  @param[in] rhs: Right hand side of multiply
 *  @returns New GaussPoly
 */
template <int D> GaussPoly<D> GaussFunc<D>::mult(const GaussFunc<D> &rhs) {
    GaussFunc<D> &lhs = *this;
    GaussPoly<D> result;
    result.multPureGauss(lhs, rhs);
    for (int d = 0; d < D; d++) {
        double newPos = result.getPos()[d];
        Polynomial lhsPoly(newPos - lhs.getPos()[d], lhs.getPow(d));
        Polynomial rhsPoly(newPos - rhs.getPos()[d], rhs.getPow(d));
        Polynomial newPoly = lhsPoly * rhsPoly;
        result.setPoly(d, newPoly);
    }
    result.setCoef(result.getCoef() * lhs.getCoef() * rhs.getCoef());
    return result;
}

/** @brief Multiply GaussFunc by scalar
 *  @param[in] c: Scalar to multiply
 *  @returns New GaussFunc
 */
template <int D> GaussFunc<D> GaussFunc<D>::mult(double c) {
    GaussFunc<D> g = *this;
    g.coef *= c;
    return g;
}

template <int D> std::ostream &GaussFunc<D>::print(std::ostream &o) const {
    auto is_array = details::are_all_equal<D>(this->getExp());

    // If all of the values in the exponential are the same only
    // one is printed, else, all of them are printed.

    o << "Coef    : " << this->getCoef() << std::endl;
    if (!is_array) {
        o << "Exp     : ";
        for (auto &alpha : this->getExp()) o << alpha << " ";
    } else {
        o << "Exp     : " << this->getExp(0) << std::endl;
    }
    o << "Pos     : ";
    for (int i = 0; i < D; i++) o << this->getPos()[i] << " ";
    o << std::endl;
    o << "Pow     : ";
    for (int i = 0; i < D; i++) o << this->getPow()[i] << " ";
    o << std::endl;
    return o;
}

/** @brief Compute Coulomb repulsion energy between two GaussFuncs
 *  @param[in] this: Left hand GaussFunc
 *  @param[in] rhs: Right hand GaussFunc
 *  @returns Coulomb energy
 *
 *  @note Both Gaussians must be normalized to unit charge
 *  \f$ \alpha = (\beta/\pi)^{D/2} \f$ for this to be correct!
 */
template <int D> double GaussFunc<D>::calcCoulombEnergy(const GaussFunc<D> &gf) const {
    NOT_IMPLEMENTED_ABORT;
}

template <> double GaussFunc<3>::calcCoulombEnergy(const GaussFunc<3> &gf) const {

    // Checking if the elements in each exponent are constant
    if (!details::are_all_equal<3>(this->getExp()) or !details::are_all_equal<3>(gf.getExp())) NOT_IMPLEMENTED_ABORT;

    // If they are constant the 0th element are assigned a value
    // and the Coulomb Energy can be calculated
    auto p = this->getExp()[0];
    auto q = gf.getExp()[0];

    double alpha = p * q / (p + q);

    const auto &Rp = this->getPos();
    const auto &Rq = gf.getPos();

    double Rx = Rp[0] - Rq[0];
    double Ry = Rp[1] - Rq[1];
    double Rz = Rp[2] - Rq[2];

    double Rpq_2 = Rx * Rx + Ry * Ry + Rz * Rz;

    BoysFunction boys(0);

    Coord<1> boysArg{alpha * Rpq_2};
    double boysFac = boys.evalf(boysArg);

    return std::sqrt(4.0 * alpha / pi) * boysFac;
}

template class GaussFunc<1>;
template class GaussFunc<2>;
template class GaussFunc<3>;
} // namespace mrcpp
