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
 * # Polynomial (interface)
 *
 * A light-weight, *affine-mapped* univariate polynomial used throughout MRCPP.
 * Internally, a polynomial is represented in the auxiliary variable
 *
 *   \f$ q = N\,x - L \f$
 *
 * where `N` is a dilation and `L` a translation. Coefficients are stored
 * in **ascending** powers of `q`, i.e. `coefs[k]` multiplies \f$q^k\f$.
 *
 * The class supports:
 * - optional finite **bounds** (via the @ref RepresentableFunction base);
 *   values outside the bounds evaluate to 0.
 * - algebra (sum, product, scalar scale) **within the same affine map**
 *   (same `N` and `L`).
 * - analytical **derivatives**, **antiderivatives**, **inner products**
 *   and **definite integrals**.
 *
 * ## Affine operations (N, L)
 * - `setDilation`, `setTranslation` overwrite the affine map.
 * - `dilate(n)` changes the current map as `N ← N*n`.
 * - `translate(l)` applies an external x-translation by `l`, which in the
 *   internal map becomes `L ← L + N*l` so that the *external* shift is by `l`.
 *
 * ## Order vs. size
 * - `size()` returns the raw length of the coefficient vector.
 * - `getOrder()` returns the highest index whose coefficient is numerically
 *   non-zero (trims trailing ~0 entries defined by `MachineZero`).
 */

#pragma once

#include <vector>

#include <Eigen/Core>

#include "RepresentableFunction.h"

namespace mrcpp {

class Polynomial : public RepresentableFunction<1, double> {
public:
    /** @name Constructors
     *  @{
     */
    /** @brief Zero polynomial of order @p k on optional bounds [a,b]. */
    Polynomial(int k = 0, const double *a = nullptr, const double *b = nullptr);
    /** @overload */
    Polynomial(int k, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(k, a.data(), b.data()) {}
    /** @brief From coefficient vector (ascending powers in q) and optional bounds. */
    Polynomial(const Eigen::VectorXd &c, const double *a = nullptr, const double *b = nullptr);
    /** @overload */
    Polynomial(const Eigen::VectorXd &c, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(c, a.data(), b.data()) {}
    /**
     * @brief Constructs the binomial expansion of \f$(x-c)^k\f$ with optional bounds.
     *
     * Coefficients are filled using the binomial theorem; the internal map is
     * initialized to the identity (`N=1, L=0`).
     */
    Polynomial(double c, int k = 0, const double *a = nullptr, const double *b = nullptr);
    /** @overload */
    Polynomial(double c, int k, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(c, k, a.data(), b.data()) {}
    /** @brief Deep copy (including bounds and affine map). */
    Polynomial(const Polynomial &poly);
    /** @brief Deep copy assignment (including bounds and affine map). */
    Polynomial &operator=(const Polynomial &poly);
    virtual ~Polynomial() = default;
    /** @} */

    /** @name Evaluation
     *  @{
     */
    /**
     * @brief Evaluate at external coordinate \f$x\f$.
     *
     * If the polynomial has active bounds, returns `0` outside the bounded
     * interval (in x). Internally evaluates the q-series with
     * \f$q = N x - L\f$.
     */
    double evalf(double x) const;
    /** @brief Convenience overload using a @ref Coord wrapper. */
    double evalf(const Coord<1> &r) const { return evalf(r[0]); }
    /** @} */

    /** @name Bounds mapped to x
     *  @{
     */
    /** @brief Lower bound in x corresponding to the internal bound in q. */
    double getScaledLowerBound() const;
    /** @brief Upper bound in x corresponding to the internal bound in q. */
    double getScaledUpperBound() const;
    /** @} */

    /** @name Norms
     *  @{
     */
    /** @brief L2-normalize on current (finite) bounds; no-op if unbounded. */
    void normalize();
    /**
     * @brief Squared L2 norm on current bounds.
     * @return \f$\|P\|^2\f$ if bounded; `-1` if unbounded.
     */
    double calcSquareNorm();
    /** @} */

    /** @name Affine map (q = N x - L)
     *  @{
     */
    double getTranslation() const { return this->L; } ///< Current L (translation in q-map).
    double getDilation()   const { return this->N; } ///< Current N (dilation in q-map).

    void setDilation(double n)    { this->N = n; }                 ///< Overwrite N.
    void setTranslation(double l) { this->L = l; }                 ///< Overwrite L.
    void dilate(double n)         { this->N *= n; }                ///< Scale N in place.
    /**
     * @brief External x-translation by @p l.
     *
     * Adjusts the internal map as \f$L \leftarrow L + N\,l\f$ so that
     * \f$q = N(x+l) - L_\text{old} = N x - (L_\text{old}-N l)\f$.
     */
    void translate(double l)      { this->L += this->N * l; }
    /** @} */

    /** @name Coefficients and order
     *  @{
     */
    int size() const { return this->coefs.size(); } ///< Raw length of the coefficient vector (q-powers).
    /**
     * @brief Highest non-negligible power (polynomial degree).
     *
     * Scans from low to high and returns the largest index whose coefficient
     * magnitude exceeds `MachineZero`. May be smaller than `size()-1`.
     */
    int getOrder() const;
    /** @brief Replace coefficients with a single zero (reset to degree 0). */
    void clearCoefs() { this->coefs = Eigen::VectorXd::Zero(1); }
    /** @brief Zero all current coefficients (preserve vector length). */
    void setZero() { this->coefs = Eigen::VectorXd::Zero(this->coefs.size()); }
    /** @brief Overwrite the coefficient vector (ascending powers in q). */
    void setCoefs(const Eigen::VectorXd &c) { this->coefs = c; }

    /** @brief Mutable access to the coefficient vector. */
    Eigen::VectorXd &getCoefs() { return this->coefs; }
    /** @brief Const access to the coefficient vector. */
    const Eigen::VectorXd &getCoefs() const { return this->coefs; }
    /** @} */

    /** @name Calculus
     *  @{
     */
    /** @brief Returns \f$ P' \f$ (derivative w.r.t. x). */
    Polynomial calcDerivative() const;
    /** @brief Returns an antiderivative \f$ Q \f$ with \f$Q(0)=0\f$. */
    Polynomial calcAntiDerivative() const;

    /** @brief In-place derivative \f$ P \leftarrow P' \f$. */
    void calcDerivativeInPlace();
    /** @brief In-place antiderivative \f$ P \leftarrow \int P\,dx \f$, constant = 0. */
    void calcAntiDerivativeInPlace();
    /** @} */

    /** @name Integration & inner product
     *  @{
     */
    /**
     * @brief Analytic definite integral \f$\int_a^b P(x)\,dx\f$.
     *
     * - If the polynomial has internal bounds, integrates over the
     *   intersection with \f$[a,b]\f$ (if `a`/`b` are provided).
     * - If unbounded, both `a` and `b` must be provided.
     */
    double integrate(const double *a = 0, const double *b = 0) const;
    /**
     * @brief Inner product \f$\langle P,Q\rangle = \int P(x)Q(x)\,dx\f$ over P's bounds.
     *
     * Requires `*this` to be bounded. The product is formed algebraically and
     * integrated over the same interval.
     */
    double innerProduct(const Polynomial &p) const;
    /** @} */

    /** @name Algebra (same affine map required)
     *  @{
     */
    /**
     * @brief Fused add: \f$ P \leftarrow P + c\,Q \f$.
     *
     * @note Both operands must have the same `(N,L)`; this is enforced in the
     * implementation and will error out if violated.
     */
    void addInPlace(double c, const Polynomial &Q);
    /** @brief Returns \f$ R = P + c\,Q \f$ (operands unchanged). */
    Polynomial add(double c, const Polynomial &Q) const;

    /** @brief Scalar product \f$ Q = c\,P \f$. */
    Polynomial operator*(double c) const;
    /**
     * @brief Polynomial product \f$ R = P\cdot Q \f$.
     *
     * @note Requires same `(N,L)` affine map in the implementation.
     */
    Polynomial operator*(const Polynomial &Q) const;

    /** @brief Sum \f$ P+Q \f$ (convenience). */
    Polynomial operator+(const Polynomial &Q) const { return add(1.0, Q); }
    /** @brief Difference \f$ P-Q \f$ (convenience). */
    Polynomial operator-(const Polynomial &Q) const { return add(-1.0, Q); }

    /** @brief In-place scalar scale: \f$ P \leftarrow c\,P \f$. */
    Polynomial &operator*=(double c);
    /**
     * @brief In-place product: \f$ P \leftarrow P\cdot Q \f$.
     *
     * @note Requires same `(N,L)` affine map in the implementation.
     */
    Polynomial &operator*=(const Polynomial &Q);
    /** @brief In-place sum: \f$ P \leftarrow P+Q \f$. */
    Polynomial &operator+=(const Polynomial &Q);
    /** @brief In-place difference: \f$ P \leftarrow P-Q \f$. */
    Polynomial &operator-=(const Polynomial &Q);
    /** @} */

protected:
    double N;              ///< Dilation in the internal map \f$q = N x - L\f$.
    double L;              ///< Translation in the internal map \f$q = N x - L\f$.
    Eigen::VectorXd coefs; ///< Coefficients for ascending powers of \f$q\f$.
};

} // namespace mrcpp