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

#pragma once

#include <Eigen/Core>
#include <vector>

#include "MRCPP/mrcpp_declarations.h"

#include "Gaussian.h"
#include "Polynomial.h"

namespace mrcpp {

/** @class GaussPoly
 *
 * @brief Polynomial–Gaussian in D dimensions (separable form).
 *
 * Concept
 * -------
 * GaussPoly represents functions of the form
 * \f[
 *   f(\mathbf x) \;=\; c \;\prod_{d=1}^{D}\; P_d(x_d - x^{(0)}_d)\,
 *                      \exp\!\big(-\alpha_d\,(x_d - x^{(0)}_d)^2\big),
 * \f]
 * i.e. a per–dimension polynomial factor times an anisotropic Gaussian.
 * The per–axis polynomials \f$P_d\f$ are stored as owned pointers
 * `Polynomial* poly[d]`.  This class derives from @ref Gaussian to reuse
 * storage for coefficient `coef`, exponents `alpha[d]`, and center `pos[d]`.
 *
 * Key features
 * ------------
 * - Exact evaluation in 1D/ND (see evalf / evalf1D).
 * - Exact L2–norm via expansion into a sum of monomial–Gaussians
 *   (@ref asGaussExp + analytic overlaps).
 * - Algebraic utilities (scalar and poly products; some are intentionally
 *   left unimplemented to avoid accidental heavy symbolic work).
 *
 * Ownership
 * ---------
 * The `poly[d]` pointers are owned by the GaussPoly instance (deep-copied
 * in copy-operations and destroyed in the destructor).
 *
 * Relation to GaussFunc
 * ---------------------
 * A @ref GaussFunc corresponds to the special case where each `P_d(t)=t^{k_d}`
 * is a monomial.  A GaussPoly can be expanded to a sum of GaussFunc terms
 * (tensor product of monomials) with @ref asGaussExp.
 */

template <int D> class GaussPoly : public Gaussian<D> {
public:
    /** @name Constructors & Lifetime
     *  @{
     */

    /** @brief Construct an isotropic GaussPoly with optional per-axis degrees.
     *
     * @param[in] alpha  Exponent parameter (isotropic): \f$ \alpha_d \equiv \alpha \f$.
     * @param[in] coef   Global amplitude \f$ c \f$.
     * @param[in] pos    Center \f$ x^{(0)} \f$ per dimension.
     * @param[in] power  Maximum polynomial degree per dimension (order of @ref Polynomial).
     *
     * Initializes each `poly[d]` as a Polynomial of degree `power[d]`
     * (if any non-zero degree is requested), otherwise keeps it nullptr.
     * The Gaussian base class stores `(coef, alpha, pos, power)`.
     */
    GaussPoly(double alpha = 0.0, double coef = 1.0, const Coord<D> &pos = {}, const std::array<int, D> &power = {});

    /** @brief Construct an anisotropic GaussPoly (per-axis exponents).
     *
     * @param[in] alpha  Per-axis exponents \f$ \{\alpha_d\}_{d=1}^D \f$.
     * @param[in] coef   Global amplitude.
     * @param[in] pos    Center per dimension.
     * @param[in] power  Maximum polynomial degree per dimension.
     *
     * Same allocation policy for `poly[d]` as in the isotropic constructor.
     */
    GaussPoly(const std::array<double, D> &alpha,
              double coef,
              const Coord<D> &pos = {},
              const std::array<int, D> &power = {});

    /** @brief Deep-copy ctor (also clones per-axis polynomials). */
    GaussPoly(const GaussPoly<D> &gp);

    /** @brief Build GaussPoly from a @ref GaussFunc (monomial×Gaussian).
     *
     * Creates per-axis polynomials equal to the corresponding monomials,
     * i.e. `P_d(t) = t^{power[d]}`.
     */
    GaussPoly(const GaussFunc<D> &gf);

    /** @brief Disable copy-assignment (explicit semantic/ownership choice). */
    GaussPoly<D> &operator=(const GaussPoly<D> &gp) = delete;

    /** @brief Polymorphic clone (deep copy). */
    Gaussian<D> *copy() const override;

    /** @brief Destructor; releases owned Polynomial pointers. */
    ~GaussPoly();

    /** @} */

    /** @name Math & Evaluation
     *  @{
     */

    /** @brief Exact L2-norm squared \f$ \|f\|_2^2 \f$.
     *
     * Implementation:
     * 1) Expand to a sum of monomial Gaussians (@ref asGaussExp).
     * 2) Sum analytic overlaps of all pairs (Obara–Saika), see
     *    `function_utils::calc_overlap`.
     */
    double calcSquareNorm() const override;

    /** @brief Evaluate \f$ f(\mathbf x) \f$ at a point (D-D). */
    double evalf(const Coord<D> &r) const override;

    /** @brief Evaluate the 1D factor in dimension `dim` at coordinate `r`.
     *
     * The convention (consistent with other classes): the global amplitude
     * `coef` is applied only in `dim==0` so that a tensor product across
     * dimensions yields the correct global amplitude once.
     */
    double evalf1D(double r, int dim) const override;

    /** @brief Expand into a sum of @ref GaussFunc terms (tensor of monomials).
     *
     * Produces \f$ \prod_d P_d \f$ as a sum of monomials and attaches the same
     * Gaussian envelope.  This is used both for integration and algebra.
     */
    GaussExp<D> asGaussExp() const override;

    /** @brief Analytic derivative in Cartesian direction `dir`.
     *
     * @note The implementation may throw/abort if not provided for GaussPoly.
     *       (The .cpp currently marks this as NOT_IMPLEMENTED.)
     */
    GaussPoly differentiate(int dir) const override;

    /** @} */

    /** @name Algebra
     *  @{
     */

    /** @brief In-place product with another GaussPoly (same center/envelope).
     *
     * @warning Not implemented in the current source (will abort if called).
     */
    void multInPlace(const GaussPoly<D> &rhs);

    /** @brief In-place product operator (delegates to @ref multInPlace). */
    void operator*=(const GaussPoly<D> &rhs) { multInPlace(rhs); }

    /** @brief Symbolic product, returns a new GaussPoly.
     *
     * @warning Not implemented in the current source (will abort if called).
     */
    GaussPoly<D> mult(const GaussPoly<D> &rhs);

    /** @brief Multiply by scalar (returns a copy). */
    GaussPoly<D> mult(double c);

    /** @brief Operator sugar for @ref mult(const GaussPoly&). */
    GaussPoly<D> operator*(const GaussPoly<D> &rhs) { return mult(rhs); }

    /** @brief Operator sugar for @ref mult(double). */
    GaussPoly<D> operator*(double c) { return mult(c); }

    /** @} */

    /** @name Accessors (per-axis polynomials)
     *  @{
     */

    /** @brief Read-only access to coefficient vector of polynomial in dim `i`. */
    const Eigen::VectorXd &getPolyCoefs(int i) const { return poly[i]->getCoefs(); }

    /** @brief Mutable access to coefficient vector of polynomial in dim `i`. */
    Eigen::VectorXd &getPolyCoefs(int i) { return poly[i]->getCoefs(); }

    /** @brief Read-only access to polynomial object in dim `i`. */
    const Polynomial &getPoly(int i) const { return *poly[i]; }

    /** @brief Mutable access to polynomial object in dim `i`. */
    Polynomial &getPoly(int i) { return *poly[i]; }

    /** @} */

    /** @name Mutators (structure/shape)
     *  @{
     */

    /** @brief Set polynomial degree in one dimension (reallocates @ref Polynomial). */
    void setPow(int d, int pow) override;

    /** @brief Set polynomial degrees in all dimensions (reallocates). */
    void setPow(const std::array<int, D> &pow) override;

    /** @brief Replace polynomial in dimension `d` with a copy of `poly`.
     *
     * Updates the stored per-axis degree to `poly.getOrder()`.
     * Ownership remains with this GaussPoly (deep copy).
     */
    void setPoly(int d, Polynomial &poly);

    /** @} */

private:
    /** @brief Owned per-axis polynomials \f$P_d\f$ (nullptr if unused). */
    Polynomial *poly[D];

    /** @brief Helper (recursive): enumerate all monomial power combinations
     *         and collect combined coefficients (raw C-array version).
     *
     * Used by @ref asGaussExp to create the full tensor expansion.  On the
     * recursion leaf it pushes:
     *  - a newly allocated `int[D]` with the current powers, and
     *  - the corresponding scalar coefficient (product of axis coefficients,
     *    times the global amplitude).
     */
    void fillCoefPowVector(std::vector<double> &coefs, std::vector<int *> &power, int pow[D], int dir) const;

    /** @brief Helper (recursive): same as above, with std::array accumulator. */
    void fillCoefPowVector(std::vector<double> &coefs,
                           std::vector<int *> &power,
                           std::array<int, D> &pow,
                           int dir) const;

    /** @brief Pretty-print (polynomial degrees, coefficients, envelope). */
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp