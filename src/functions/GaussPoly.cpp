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

#include "GaussPoly.h"

#include <vector>

#include "function_utils.h"
#include "GaussExp.h"
#include "GaussFunc.h"
#include "utils/Printer.h"
#include "utils/details.h"

using namespace Eigen;

namespace mrcpp {

/** @returns New GaussPoly object
 *  @param[in] beta: Exponent, \f$ e^{-\beta r^2} \f$
 *  @param[in] alpha: Coefficient, \f$ \alpha e^{-r^2} \f$
 *  @param[in] pos: Position \f$ (x - pos[0]), (y - pos[1]), ... \f$
 *  @param[in] pow: Max polynomial degree, \f$ P_0(x), P_1(y), ... \f$
 *
 *  High-level:
 *  -----------
 *  GaussPoly<D> represents a separable polynomial-times-Gaussian:
 *      f(x) = coef * Π_d Poly_d(x_d - pos[d]) * exp( -alpha[d] (x_d - pos[d])^2 ).
 *  The per-axis polynomials are stored as pointers (Polynomial* poly[d]).
 *  Here we allocate those polynomials (if a non-zero power is requested),
 *  using the degree from `power[d]` passed to the Gaussian base ctor.
 */
template <int D>
GaussPoly<D>::GaussPoly(double beta, double alpha, const Coord<D> &pos, const std::array<int, D> &power)
        : Gaussian<D>(beta, alpha, pos, power) {
    for (auto d = 0; d < D; d++) {
        // If overall 'power' array is not the all-zero sentinel, create a poly
        // of the requested degree for this axis. Otherwise leave pointer null.
        if (power != std::array<int, D>{}) {
            this->poly[d] = new Polynomial(this->power[d]);
        } else {
            this->poly[d] = nullptr;
        }
    }
}

/** @brief Anisotropic exponent ctor (per-axis beta).
 *  Same allocation logic for the per-axis polynomials as above.
 */
template <int D>
GaussPoly<D>::GaussPoly(const std::array<double, D> &beta,
                        double alpha,
                        const Coord<D> &pos,
                        const std::array<int, D> &pow)
        : Gaussian<D>(beta, alpha, pos, pow) {
    for (auto d = 0; d < D; d++) {
        if (pow != std::array<int, D>{}) {
            this->poly[d] = new Polynomial(this->power[d]);
        } else {
            this->poly[d] = nullptr;
        }
    }
}

/** @brief Copy-construct with deep copies of the per-axis polynomials. */
template <int D>
GaussPoly<D>::GaussPoly(const GaussPoly<D> &gp)
        : Gaussian<D>(gp) {
    for (int d = 0; d < D; d++) { poly[d] = new Polynomial(gp.getPoly(d)); }
}

/** @brief Construct a GaussPoly from a GaussFunc (monomial×Gaussian).
 *
 *  Effect:
 *  -------
 *  For each axis d, we create a polynomial of degree equal to the monomial power
 *  in that dimension, and set it to the *monomial basis* e_d(t) = t^{power[d]}
 *  (i.e., coefficient vector with a single 1 at index = power[d]).
 */
template <int D>
GaussPoly<D>::GaussPoly(const GaussFunc<D> &gf)
        : Gaussian<D>(gf) {
    for (int d = 0; d < D; d++) {
        int order = this->getPower(d);
        poly[d] = new Polynomial(order);
        VectorXd coefs = VectorXd::Zero(order + 1);
        coefs[order] = 1.0;               // t^{order}
        poly[d]->setCoefs(coefs);
        // poly[d]->unsetBounds();
    }
}

/** @brief Delete owned Polynomial objects. */
template <int D> GaussPoly<D>::~GaussPoly() {
    for (int i = 0; i < D; i++) { delete poly[i]; }
}

/** @brief Virtual clone (deep copy). */
template <int D> Gaussian<D> *GaussPoly<D>::copy() const {
    auto *gauss = new GaussPoly<D>(*this);
    return gauss;
}

/** @brief Exact L2 norm squared by expanding to GaussExp and summing overlaps.
 *
 *  Algorithm:
 *  ----------
 *   1) Expand this GaussPoly into a sum of GaussFunc terms (asGaussExp()).
 *   2) Sum ⟨g_i | g_j⟩ over all pairs using the analytic overlap routine
 *      function_utils::calc_overlap (Obara–Saika recurrences).
 */
template<int D> double GaussPoly<D>::calcSquareNorm() const {
    GaussExp<D> this_exp = this->asGaussExp();
    double norm = 0.0;
    for (int i = 0; i < this_exp.size(); i++) {
        auto func_i = static_cast<GaussFunc<D> &>(this_exp.getFunc(i));
        for (int j = 0; j < this_exp.size(); j++) {
            auto func_j = static_cast<GaussFunc<D> &>(this_exp.getFunc(j));
            norm += function_utils::calc_overlap(func_i, func_j);
        }
    }
    return norm;
}

/** @brief Evaluate f(r) = coef * Π_d poly_d(r_d - pos[d]) * exp(-Σ_d alpha[d](r_d-pos[d])^2)
 *
 *  Notes:
 *  ------
 *  - Optional *screening*: if enabled, points outside the [A,B] box give 0.
 *  - The polynomial is evaluated in *shifted* coordinate q = r_d - pos[d].
 */
template <int D> double GaussPoly<D>::evalf(const Coord<D> &r) const {
    if (this->getScreen()) {
        for (int d = 0; d < D; d++) {
            if (r[d] < this->A[d] or r[d] > this->B[d]) { return 0.0; }
        }
    }
    double q2 = 0.0, p2 = 1.0;
    for (int d = 0; d < D; d++) {
        // assert(this->poly[d]->getCheckBounds() == false);
        double q = r[d] - this->pos[d];
        q2 += this->alpha[d] * q * q;           // accumulate quadratic exponent
        p2 *= poly[d]->evalf(r[d] - this->pos[d]); // polynomial factor in dim d
    }
    return this->coef * p2 * std::exp(-q2);
}

/** @brief Evaluate the *1D* factor in dimension d at coordinate r.
 *
 *  Implementation detail:
 *  ----------------------
 *  For efficiency, only dimension d=0 gets the *full* global coefficient.
 *  Other dimensions return the pure 1D factor with amplitude 1.0. This is a
 *  deliberate convention to avoid taking the d-th root of the coefficient when
 *  forming tensor products; callers multiply across dims and obtain the correct
 *  full amplitude once (from d==0).
 */
template <int D> double GaussPoly<D>::evalf1D(const double r, int d) const {
    // NOTE!
    //     This function evaluation will give the first dimension the full coef
    //     amplitude, leaving all other directions with amplitude 1.0. This is to
    //     avoid expensive d-root evaluation when distributing the amplitude
    //     equally to all dimensions.

    if (this->getScreen()) {
        if ((r < this->A[d]) or (r > this->B[d])) { return 0.0; }
    }
    // assert(this->poly[d]->getCheckBounds() == false);
    double q2 = 0.0, p2 = 1.0;
    double q = (r - this->pos[d]);
    q2 += q * q;
    p2 *= poly[d]->evalf(q);
    if (d == 0) { p2 *= this->coef; }      // apply global amplitude once
    return p2 * std::exp(-this->alpha[d] * q2);
}

/** @brief Expand a polynomial×Gaussian into a sum of pure Gaussians (GaussExp).
 *
 *  Idea:
 *  -----
 *  Each per-axis polynomial Poly_d(t) = Σ_{k=0}^{p_d} c_{d,k} t^k can be viewed
 *  as a linear combination of *monomial* GaussFuncs: (x-pos[d])^k * exp(-α_d t^2).
 *  The D-dimensional product of polynomials expands into a tensor product of
 *  monomials across dimensions. This routine enumerates all combinations of
 *  powers (k_0,...,k_{D-1}), multiplies coefficients Π_d c_{d,k_d}, and emits
 *  corresponding GaussFunc terms into a GaussExp.
 *
 *  Implementation:
 *  ---------------
 *  - nTerms = Π_d (power[d] + 1).
 *  - fillCoefPowVector(...) recursively builds:
 *      * `coefs[i]` = Π_d c_{d, pow_d(i)} * global coef
 *      * `power[i]` = array<int,D> of the per-axis monomial powers
 *  - For each nonzero coefficient, create GaussFunc(alpha, coef, pos, pow).
 */
template <int D> GaussExp<D> GaussPoly<D>::asGaussExp() const {
    std::array<int, D> pow;
    std::array<double, D> pos;
    auto alpha = this->getExp();

    int nTerms = 1;
    for (int d = 0; d < D; d++) {
        nTerms *= (this->getPower(d) + 1);
        pos[d] = this->getPos()[d];
    }

    std::vector<double> coefs;
    std::vector<int *> power;

    fillCoefPowVector(coefs, power, pow, D);

    GaussExp<D> gexp;
    for (int i = 0; i < nTerms; i++) {
        double coef = coefs[i];
        for (int d = 0; d < D; d++) pow[d] = power[i][d];
        if (coef != 0.0) {
            GaussFunc<D> gFunc(alpha, coef, pos, pow);
            gexp.append(gFunc);
        }
    }
    for (auto &i : power) { delete[] i; }
    return gexp;
}

/** @brief Analytic derivative (not implemented for GaussPoly). */
template <int D> GaussPoly<D> GaussPoly<D>::differentiate(int dir) const {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief In-place product with another GaussPoly (not implemented). */
template <int D> void GaussPoly<D>::multInPlace(const GaussPoly<D> &rhs) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief Recursive helper: enumerate all power combinations; collect coefficients.
 *
 *  Version 1: temporary raw `int pow[D]` buffer.
 *
 *  On the recursion leaf (dir==0 processed), allocate a new array `newPow[d]`
 *  storing the tuple of powers; compute the scalar coefficient as:
 *      coef = global_coef * Π_d poly_d->coefs[ pow[d] ]
 *  and push both into the output vectors.
 */
template <int D>
void GaussPoly<D>::fillCoefPowVector(std::vector<double> &coefs, std::vector<int *> &power, int pow[D], int dir) const {
    dir--;
    for (int i = 0; i < this->getPower(dir) + 1; i++) {
        pow[dir] = i;
        if (dir > 0) {
            fillCoefPowVector(coefs, power, pow, dir);
        } else {
            auto *newPow = new int[D];
            double coef = 1.0;
            for (int d = 0; d < D; d++) {
                newPow[d] = pow[d];
                coef *= this->getPolyCoefs(d)[pow[d]];
            }
            coef *= this->getCoef();
            power.push_back(newPow);
            coefs.push_back(coef);
        }
    }
}

/** @brief Recursive helper: same as above, with std::array<int,D> accumulator. */
template <int D>
void GaussPoly<D>::fillCoefPowVector(std::vector<double> &coefs,
                                     std::vector<int *> &power,
                                     std::array<int, D> &pow,
                                     int dir) const {
    dir--;
    for (int i = 0; i < this->getPower(dir) + 1; i++) {
        pow[dir] = i;
        if (dir > 0) {
            fillCoefPowVector(coefs, power, pow, dir);
        } else {
            auto *newPow = new int[D];
            double coef = 1.0;
            for (int d = 0; d < D; d++) {
                newPow[d] = pow[d];
                coef *= this->getPolyCoefs(d)[pow[d]];
            }
            coef *= this->getCoef();
            power.push_back(newPow);
            coefs.push_back(coef);
        }
    }
}

/** @brief Product of two GaussPoly (symbolic) — currently not implemented. */
template <int D> GaussPoly<D> GaussPoly<D>::mult(const GaussPoly<D> &rhs) {
    NOT_IMPLEMENTED_ABORT;
    /*
    GaussPoly<D> &lhs = *this;
    GaussPoly<D> result;
    result.multPureGauss(lhs, rhs);
    for (int d = 0; d < D; d++) {
        double newPos = result.getPos()[d];
        int lhsPow = lhs.getPower(d);
        Polynomial lhsPoly(lhsPow);
        lhsPoly.clearCoefs();
        for (int p = 0; p <= lhsPow; p++) {
            Polynomial tmpPoly(newPos - lhs.getPos()[d], p);
            tmpPoly *= lhs.getPolyCoefs(d)[p];
            lhsPoly += tmpPoly;
        }

        int rhsPow = rhs.getPower(d);
        Polynomial rhsPoly(rhsPow);
        rhsPoly.clearCoefs();
        for (int p = 0; p <= rhsPow; p++) {
            Polynomial tmpPoly(newPos - rhs.getPos()[d], p);
            tmpPoly *= rhs.getPolyCoefs(d)[p];
            rhsPoly += tmpPoly;
        }
        Polynomial newPoly = lhsPoly * rhsPoly;
        result.setPoly(d, newPoly);
    }
    result.setCoef(result.getCoef() * lhs.getCoef() * rhs.getCoef());
    return result;
    */
}

/** @brief Multiply by a scalar (returns a copy). */
template <int D> GaussPoly<D> GaussPoly<D>::mult(double c) {
    GaussPoly<D> g = *this;
    g.coef *= c;
    return g;
}

/** @brief Set polynomial degree *and* allocate a new Polynomial in dim d. */
template <int D> void GaussPoly<D>::setPow(int d, int pow) {
    if (poly[d] != nullptr) { delete poly[d]; }
    poly[d] = new Polynomial(pow);
}

/** @brief Set polynomial degrees in all dims and allocate new polynomials. */
template <int D> void GaussPoly<D>::setPow(const std::array<int, D> &pow) {
    for (int d = 0; d < D; d++) {
        if (poly[d] != nullptr) { delete poly[d]; }
        poly[d] = new Polynomial(pow[d]);
    }
}

/** @brief Replace the polynomial in a given dimension and update degree.
 *
 *  Ownership:
 *  ----------
 *  This class owns its per-axis Polynomial pointers. We take a *copy* of the
 *  passed polynomial to keep ownership consistent and update power[d] to match
 *  the new polynomial order.
 */
template <int D> void GaussPoly<D>::setPoly(int d, Polynomial &poly) {
    if (this->poly[d] != nullptr) { delete this->poly[d]; }
    this->poly[d] = new Polynomial(poly);
    this->power[d] = poly.getOrder();
}

/** @brief Pretty-print parameters, including per-axis polynomial coefficients. */
template <int D> std::ostream &GaussPoly<D>::print(std::ostream &o) const {
    auto is_array = details::are_all_equal<D>(this->getExp());

    // If all exponents are identical, print a single value; else print the array.
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
    for (int i = 0; i < D; i++) o << "Poly[" << i << "] : " << this->getPolyCoefs(i).transpose() << std::endl;
    return o;
}

// Explicit template instantiations
template class GaussPoly<1>;
template class GaussPoly<2>;
template class GaussPoly<3>;

} // namespace mrcpp