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
 * @file GaussFunc.cpp
 *
 * @brief Implementation of @c GaussFunc<D>, a single Cartesian Gaussian (possibly
 *        multiplied by a coordinate power) in D dimensions.
 *
 * Model
 * -----
 * A term has the form
 *   f(r) = c * Π_{d=0}^{D-1} (x_d - R_d)^{p_d} * exp( -α_d (x_d - R_d)^2 ),
 * with scalar coefficient c, center R, exponents α (per axis), and integer powers p
 * (Cartesian angular momenta). Many operations here are separable over dimensions.
 *
 * Highlights
 * ----------
 * - @ref evalf computes the value with optional screening (box truncation).
 * - @ref calcSquareNorm uses 1D closed forms and multiplies across axes.
 * - @ref differentiate returns a @ref GaussPoly<D> (Gaussian times polynomial),
 *   using the product rule on (x-R)^p * exp(-α (x-R)^2).
 * - @ref mult multiplies two @c GaussFunc into a @c GaussPoly by “completing the square”
 *   (handled by @c GaussPoly::multPureGauss) and then combining the two polynomials
 *   created by shifting to the new center.
 * - @ref calcCoulombEnergy (D=3 specialization) uses Boys F_0 and assumes isotropic
 *   exponents for both Gaussians.
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

/**
 * @brief Polymorphic deep copy (virtual constructor idiom).
 * @return Newly allocated @c GaussFunc with identical parameters.
 */
template <int D> Gaussian<D> *GaussFunc<D>::copy() const {
    auto *gauss = new GaussFunc<D>(*this);
    return gauss;
}

/**
 * @brief Pointwise evaluation of the Gaussian (with optional polynomial factor).
 *
 * Steps
 * -----
 * 1) If screening is enabled, immediately return 0 if any coordinate lies outside
 *    the precomputed box [A[d], B[d]].
 * 2) Accumulate q2 = Σ α_d (x_d - R_d)^2 (the exponent argument),
 *    and p2 = Π (x_d - R_d)^{p_d} (the Cartesian polynomial).
 * 3) Return c * p2 * exp(-q2).
 */
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

/**
 * @brief 1D evaluation of the d-th component only (factor for separable product).
 *
 * This returns (x - R_d)^{p_d} * exp(-α_d (x - R_d)^2) times @c coef if d==0
 * (the overall scalar is stored redundantly only on one axis when factoring).
 * Screening is applied on that axis if enabled.
 */
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

/**
 * @brief Squared L2 norm  ||f||^2 = ∫ |f|^2 d r  (separable product of 1D integrals).
 *
 * For one axis (drop subscript d for brevity):
 *   ∫ (x-R)^{2p} exp(-2α (x-R)^2) dx
 * = sqrt(pi / (2α)) * [(2p-1)!!] / (2α)^p,
 * which is implemented via a simple descending product:
 *    sq_norm = Π_{i odd from (2p-1) down to 1} i / (2α).
 * The D-dimensional norm is the product over axes, multiplied by coef^2.
 */
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

/**
 * @brief Convert a single @c GaussFunc into a length-1 @c GaussExp.
 *
 * Useful when operations expect an expansion (e.g., in norm cross-terms).
 */
template<int D> GaussExp<D> GaussFunc<D>::asGaussExp() const {
    GaussExp<D> gexp;
    gexp.append(*this);
    return gexp;
}

/**
 * @brief Derivative along axis @p dir, returning a @c GaussPoly<D>.
 *
 * In 1D:
 *   d/dx [(x-R)^p e^{-α(x-R)^2}] =
 *     p (x-R)^{p-1} e^{-α...}  + (x-R)^p * (-2α)(x-R) e^{-α...}
 *   = [ p (x-R)^{p-1}  - 2α (x-R)^{p+1} ] e^{-α...}
 *
 * We therefore create a new polynomial of degree (p+1) with two nonzero
 * coefficients at (p-1) and (p+1). Other axes carry over unchanged.
 */
template <int D> GaussPoly<D> GaussFunc<D>::differentiate(int dir) const {
    GaussPoly<D> result(*this);
    int oldPow = this->getPower(dir);

    Polynomial newPoly(oldPow + 1);
    newPoly.getCoefs()[oldPow + 1] = -2.0 * this->getExp()[dir];
    if (oldPow > 0) { newPoly.getCoefs()[oldPow - 1] = oldPow; }
    result.setPoly(dir, newPoly);
    ;
    return result;
}

/**
 * @brief In-place multiplication by another @c GaussFunc with the SAME center.
 *
 * Preconditions
 * -------------
 * - The two Gaussians must share identical centers in every axis.
 *
 * Effect
 * ------
 * - Exponents add: α_new = α_lhs + α_rhs.
 * - Powers add:    p_new = p_lhs + p_rhs.
 * - Coefficients multiply: c_new = c_lhs * c_rhs.
 *
 * This keeps the center unchanged and avoids creating polynomials.
 */
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
    for (int d = 0; d < D; d++) { newPow[d] = lhs.getPower(d) + rhs.getPower(d); }
    this->setCoef(newCoef);
    this->setExp(newExp);
    this->setPow(newPow);
}

/** @brief Multiply two GaussFuncs
 *  @param[in] this: Left hand side of multiply
 *  @param[in] rhs: Right hand side of multiply
 *  @returns New GaussPoly
 *
 * Algorithm
 * ---------
 * 1) “Complete the square”: the product of two Gaussians is a (shifted) Gaussian
 *    with combined exponent and a new center (weighted by exponents). This part is
 *    delegated to @c GaussPoly::multPureGauss, which sets the new Gaussian envelope
 *    (position, exponents, and a prefactor).
 * 2) Each original polynomial factor (x-R)^p is re-expressed relative to the new
 *    center R_new: (x-R) = (x-R_new) + (R_new - R). We therefore have two polynomials
 *    per axis; they are multiplied to obtain the combined polynomial for that axis.
 * 3) Multiply in the original scalar coefficients c_lhs * c_rhs.
 */
template <int D> GaussPoly<D> GaussFunc<D>::mult(const GaussFunc<D> &rhs) {
    GaussFunc<D> &lhs = *this;
    GaussPoly<D> result;
    result.multPureGauss(lhs, rhs);
    for (int d = 0; d < D; d++) {
        double newPos = result.getPos()[d];
        Polynomial lhsPoly(newPos - lhs.getPos()[d], lhs.getPower(d));
        Polynomial rhsPoly(newPos - rhs.getPos()[d], rhs.getPower(d));
        Polynomial newPoly = lhsPoly * rhsPoly;
        result.setPoly(d, newPoly);
    }
    result.setCoef(result.getCoef() * lhs.getCoef() * rhs.getCoef());
    return result;
}

/** @brief Multiply GaussFunc by scalar (returns a copy with scaled coefficient). */
template <int D> GaussFunc<D> GaussFunc<D>::mult(double c) {
    GaussFunc<D> g = *this;
    g.coef *= c;
    return g;
}

/**
 * @brief Pretty-printer for a Gaussian term.
 *
 * Prints:
 *  - Coef
 *  - Exp: either a single value if all α_d are equal, or all components.
 *  - Pos: center coordinates
 *  - Pow: integer powers per axis
 */
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
    for (int i = 0; i < D; i++) o << this->getPower()[i] << " ";
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
 *
 *  General D is not implemented here; see the D=3 specialization below.
 */
template <int D> double GaussFunc<D>::calcCoulombEnergy(const GaussFunc<D> &gf) const {
    NOT_IMPLEMENTED_ABORT;
}

/**
 * @brief Coulomb energy for 3D isotropic Gaussians using Boys F_0.
 *
 * Preconditions
 * -------------
 * - Both Gaussians must have isotropic exponents (α_x = α_y = α_z), verified via
 *   @c details::are_all_equal<3>.
 *
 * Formula
 * -------
 * With exponents p and q, α = p q / (p + q), separation R = |R_p - R_q|,
 * the Coulomb interaction is:
 *   E = sqrt( 4 α / π ) * F_0( α R^2 ),
 * where F_0 is the order-zero Boys function. The code constructs a @c BoysFunction(0)
 * and evaluates it at α R^2.
 */
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

// Explicit template instantiations
template class GaussFunc<1>;
template class GaussFunc<2>;
template class GaussFunc<3>;
} // namespace mrcpp