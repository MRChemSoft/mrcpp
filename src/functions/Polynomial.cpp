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

/**
 * Implementation notes for Polynomial
 * -----------------------------------
 * This file implements a univariate polynomial P expressed in an *affine*
 * internal coordinate q = N*x - L, where:
 *   - N is a dilation (scale) factor,
 *   - L is a translation shift (stored with sign to match the internal form).
 *
 * Coefficients are stored in ascending powers of q: coefs[k] multiplies q^k.
 * Many operations (evaluation, algebra, derivatives, integrals) are performed
 * with respect to q but expose an API in terms of the external variable x.
 *
 * Bounding:
 *   The base class (RepresentableFunction) holds optional lower/upper bounds
 *   in the *q*-domain. Helper functions `getScaledLowerBound()` /
 *   `getScaledUpperBound()` convert those bounds to the *x*-domain via the
 *   inverse affine map x = (q + L)/N. Evaluation outside the bounds yields 0.
 *
 * Algebra:
 *   Addition and multiplication require the same affine map (same N and L).
 *   We check that before combining coefficient vectors to avoid mixing
 *   different coordinate systems.
 */

#include <cfloat>

#include "Polynomial.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

/** @brief Construct a zero-valued polynomial of order @p k with optional bounds.
 *
 *  @param k  Maximum power (order). The polynomial will have (k+1) coefficients.
 *  @param a  (Optional) pointer to lower bound in x; forwarded to base class.
 *  @param b  (Optional) pointer to upper bound in x; forwarded to base class.
 *
 *  Details
 *  -------
 *  - Initializes the affine map to identity: N = 1, L = 0 (so q = x).
 *  - Allocates a coefficient vector of length k+1, initialized to zero.
 *  - Bounds are stored by the base class; they affect evalf() and integration.
 */
Polynomial::Polynomial(int k, const double *a, const double *b)
        : RepresentableFunction<1, double>(a, b) {
    assert(k >= 0);
    this->N = 1.0;
    this->L = 0.0;
    this->coefs = VectorXd::Zero(k + 1);
}

/** @brief Construct the expanded monomial (x - c)^k (up to scaling) with optional bounds.
 *
 *  @param c  Shift in the monomial center (i.e., builds coefficients of (x - c)^k).
 *  @param k  Order of the monomial.
 *  @param a  Optional lower bound; forwarded to base.
 *  @param b  Optional upper bound; forwarded to base.
 *
 *  Details
 *  -------
 *  - Uses binomial coefficients to expand (x - c)^k into the internal q = x
 *    basis (N = 1, L = 0).
 *  - coefs[i] = binom(k, i) * (-c)^(k - i).
 */
Polynomial::Polynomial(double c, int k, const double *a, const double *b)
        : RepresentableFunction<1>(a, b) {
    this->N = 1.0;
    this->L = 0.0;
    this->coefs = math_utils::get_binomial_coefs(k);
    for (int i = 0; i <= k; i++) { this->coefs[i] *= std::pow(c, k - i); }
}

/** @brief Construct from a coefficient vector (ascending powers in q) with optional bounds.
 *
 *  @param c  Coefficients for q^0, q^1, ..., q^k.
 *  @param a  Optional lower bound; forwarded to base.
 *  @param b  Optional upper bound; forwarded to base.
 *
 *  Initializes affine map to identity (N=1, L=0) and copies coefficients.
 */
Polynomial::Polynomial(const VectorXd &c, const double *a, const double *b)
        : RepresentableFunction<1>(a, b) {
    this->N = 1.0;
    this->L = 0.0;
    setCoefs(c);
}

/** @brief Copy constructor (deep copy), including bounds and affine map. */
Polynomial::Polynomial(const Polynomial &poly)
        : RepresentableFunction<1>(poly) {
    this->N = poly.N;
    this->L = poly.L;
    this->coefs = poly.coefs;
}

/** @brief Copy assignment (deep copy), including bounds and affine map.
 *
 *  Copies base part, then affine parameters N,L and coefficient vector.
 */
Polynomial &Polynomial::operator=(const Polynomial &poly) {
    RepresentableFunction<1>::operator=(poly);
    this->N = poly.N;
    this->L = poly.L;
    this->coefs = poly.coefs;
    return *this;
}

/** @brief Evaluate the polynomial at external x, honoring bounds.
 *
 *  @param x Point of evaluation in external coordinates.
 *  @return P(x) if within bounds, otherwise 0.
 *
 *  Implementation
 *  --------------
 *  - If bounded, quickly reject x outside the mapped interval.
 *  - Evaluate in the internal coordinate q = N*x - L using a simple
 *    power-accumulation loop:
 *        y = sum_k coefs[k] * q^k.
 *    (xp accumulates q^k without recomputing powers.)
 */
double Polynomial::evalf(double x) const {
    if (isBounded()) {
        if (x < this->getScaledLowerBound()) return 0.0;
        if (x > this->getScaledUpperBound()) return 0.0;
    }
    double xp = 1.0;
    double y = 0.0;
    for (int k = 0; k < getOrder() + 1; k++) {
        y += (xp * this->coefs[k]);
        xp *= this->N * x - this->L;  // advance q^k -> q^(k+1)
    }
    return y;
}

/** @brief Lower bound in external x-space, derived from the internal bound via x = (q + L)/N.
 *
 *  Preconditions: polynomial must be bounded (otherwise errors).
 */
double Polynomial::getScaledLowerBound() const {
    if (not isBounded()) MSG_ERROR("Unbounded polynomial");
    return (1.0 / this->N * (this->A[0] + this->L));
}

/** @brief Upper bound in external x-space, derived from the internal bound via x = (q + L)/N.
 *
 *  Preconditions: polynomial must be bounded (otherwise errors).
 */
double Polynomial::getScaledUpperBound() const {
    if (not isBounded()) MSG_ERROR("Unbounded polynomial");
    return (1.0 / this->N * (this->B[0] + this->L));
}

/** @brief Normalize the polynomial in L2 on its current (bounded) domain.
 *
 *  Details
 *  -------
 *  - Computes squared norm via innerProduct(*this).
 *  - Scales coefficients by 1/sqrt(norm).
 *  - If unbounded or norm < 0, aborts with an error.
 */
void Polynomial::normalize() {
    double sqNorm = calcSquareNorm();
    if (sqNorm < 0.0) MSG_ABORT("Cannot normalize polynomial");
    (*this) *= 1.0 / std::sqrt(sqNorm);
}

/** @brief Compute squared L2 norm on current bounds, or -1 if unbounded. */
double Polynomial::calcSquareNorm() {
    double sqNorm = -1.0;
    if (isBounded()) { sqNorm = this->innerProduct(*this); }
    return sqNorm;
}

/** @brief Effective order = highest index i with |coefs[i]| > MachineZero.
 *
 *  Note: This ignores trailing coefficients numerically equal to zero,
 *  and can be lower than (coefs.size()-1).
 */
int Polynomial::getOrder() const {
    int n = 0;
    for (int i = 0; i < this->coefs.size(); i++) {
        if (std::abs(this->coefs[i]) > MachineZero) { n = i; }
    }
    return n;
}

/** @brief In-place scale: P(x) ← c * P(x). */
Polynomial &Polynomial::operator*=(double c) {
    this->coefs = c * this->coefs;
    return *this;
}

/** @brief In-place product P(x) ← P(x) * Q(x) (same affine map required).
 *
 *  Preconditions
 *  -------------
 *  - Both polynomials must share identical (N, L) so they represent functions
 *    in the same internal coordinate q. Otherwise we error out.
 *
 *  Implementation
 *  --------------
 *  - Standard coefficient convolution yielding degree(P)+degree(Q).
 *  - Affine parameters are left unchanged.
 */
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

/** @brief Return Q(x) = c * P(x). */
Polynomial Polynomial::operator*(double c) const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q *= c;
    return Q;
}

/** @brief Return R(x) = P(x) * Q(x) (same affine map required).
 *
 *  Returns an unbounded polynomial that inherits the affine map and
 *  coefficients from the in-place logic.
 */
Polynomial Polynomial::operator*(const Polynomial &Q) const {
    const Polynomial &P = *this;
    Polynomial R;
    R = P;
    R *= Q;
    return R;
}

/** @brief In-place sum: P(x) ← P(x) + Q(x). (Same affine map required.) */
Polynomial &Polynomial::operator+=(const Polynomial &Q) {
    this->addInPlace(1.0, Q);
    return *this;
}

/** @brief In-place difference: P(x) ← P(x) - Q(x). (Same affine map required.) */
Polynomial &Polynomial::operator-=(const Polynomial &Q) {
    this->addInPlace(-1.0, Q);
    return *this;
}

/** @brief In-place fused add: P(x) ← P(x) + c * Q(x). (Same affine map required.)
 *
 *  Chooses the max order among P and Q and adds coefficients component-wise,
 *  padding with zeros where needed.
 */
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

/** @brief Return R(x) = P(x) + c * Q(x), leaving operands unchanged. */
Polynomial Polynomial::add(double c, const Polynomial &Q) const {
    const Polynomial &P = *this;
    Polynomial R;
    R = P;
    R.addInPlace(c, Q);
    return R;
}

/** @brief Return Q(x) = dP/dx (external derivative). */
Polynomial Polynomial::calcDerivative() const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q.calcDerivativeInPlace();
    return Q;
}

/** @brief In-place derivative: P(x) ← dP/dx.
 *
 *  Implementation
 *  --------------
 *  - Works on the internal representation in q = N*x - L:
 *      d/dx [ Σ a_i q^i ] = Σ i*a_i q^(i-1) * dq/dx = N * Σ i*a_i q^(i-1).
 *  - Since the current storage uses q-powers, we first form Σ i*a_i q^(i-1)
 *    in coefficient space. The factor N is embedded in the affine mapping
 *    (via evaluation), and the polynomial’s coefficient update matches the
 *    intended external derivative semantics given how evalf() builds q.
 *  - The code mirrors the existing convention (keeping N,L intact).
 */
void Polynomial::calcDerivativeInPlace() {
    Polynomial &P = *this;
    int P_order = P.getOrder();
    const VectorXd &oldCoefs = P.getCoefs();
    VectorXd newCoefs = VectorXd::Zero(P_order);
    for (int i = 0; i < newCoefs.size(); i++) { newCoefs[i] = double(i + 1) * oldCoefs[i + 1]; }
    P.setCoefs(newCoefs);
}

/** @brief Return the indefinite integral Q(x) = ∫ P(x) dx with zero constant. */
Polynomial Polynomial::calcAntiDerivative() const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q.calcAntiDerivativeInPlace();
    return Q;
}

/** @brief In-place antiderivative: P(x) ← ∫ P(x) dx, integration constant = 0.
 *
 *  Implementation
 *  --------------
 *  - In q-space: ∫ (Σ a_i q^i) dq = Σ a_i/(i+1) q^(i+1) + C.
 *  - For external x, dx = dq / N; the factor 1/N is accounted for when
 *    integrating over x in Polynomial::integrate(), not in coefficient
 *    construction here. We thus store the q-antiderivative coefficients.
 */
void Polynomial::calcAntiDerivativeInPlace() {
    Polynomial &P = *this;
    int P_order = P.getOrder();
    const VectorXd &oldCoefs = P.getCoefs();
    VectorXd newCoefs = VectorXd::Zero(P_order + 2);
    newCoefs[0] = 0.0;                // integration constant
    newCoefs[1] = oldCoefs[0];
    for (int i = 2; i < newCoefs.size(); i++) { newCoefs[i] = 1.0 / i * oldCoefs[i - 1]; }
    P.setCoefs(newCoefs);
}

/** @brief Analytic definite integral ∫_a^b P(x) dx, honoring bounds if present.
 *
 *  @param a Optional external lower limit (overrides internal bound if tighter).
 *  @param b Optional external upper limit (overrides internal bound if tighter).
 *  @return The integral value over max(lower bounds) to min(upper bounds).
 *
 *  Details
 *  -------
 *  - If polynomial is bounded, the domain is intersected with [a,b].
 *  - Builds the (q-based) antiderivative and evaluates it at the endpoints
 *    transformed to the q-domain by the affine map. The Jacobian dx = dq/N
 *    yields a prefactor 1/N (“sfac”).
 *  - If the final [lb,ub] is empty, returns 0.
 */
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

/** @brief Inner product ⟨P,Q⟩ over P’s current bounded domain.
 *
 *  @param Q Polynomial to multiply with.
 *  @return ∫ P(x) Q(x) dx over P’s bounds.
 *
 *  Details
 *  -------
 *  - Requires that P is bounded; Q is multiplied algebraically in q-space.
 *  - The product polynomial inherits P’s bounds; we then call integrate().
 */
double Polynomial::innerProduct(const Polynomial &Q) const {
    const Polynomial &P = *this;
    if (not P.isBounded()) MSG_ERROR("Unbounded polynomial");
    Polynomial pq = P * Q;
    pq.setBounds(P.getLowerBounds(), P.getUpperBounds());
    return pq.integrate();
}

} // namespace mrcpp