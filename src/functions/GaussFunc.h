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

#pragma once

#include <Eigen/Core>

#include "Gaussian.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/** @class GaussFunc
 *  @tparam D Spatial dimension (1,2,3,…).
 *
 *  @brief Single Cartesian Gaussian primitive (optionally with monomial powers)
 *  in D dimensions.
 *
 *  Mathematical form
 *  -----------------
 *  In D dimensions the function is separable:
 *  \f[
 *     G(\mathbf{x})
 *     = \alpha \prod_{d=0}^{D-1} (x_d - R_d)^{p_d}\,
 *       \exp\!\big(-\beta_d\,(x_d - R_d)^2\big),
 *  \f]
 *  where:
 *  - \f$ \alpha \f$ is a scalar coefficient (amplitude),
 *  - \f$ \mathbf{R} = (R_0,\dots,R_{D-1}) \f$ is the center,
 *  - \f$ \mathbf{p} = (p_0,\dots,p_{D-1}) \f$ are non-negative integers (Cartesian powers),
 *  - \f$ \boldsymbol{\beta} = (\beta_0,\dots,\beta_{D-1}) \f$ are positive exponents; they
 *    can be isotropic (\f$\beta_d=\beta\f$) or anisotropic (per-axis).
 *
 *  Relationship to @ref Gaussian
 *  -----------------------------
 *  This class *derives* from @ref Gaussian<D>, which stores the common state
 *  (coefficient, center, exponents, powers) and provides a polymorphic interface.
 *  @c GaussFunc implements operations specific to “pure Gaussian × monomial”
 *  terms (e.g., evaluation, in-place multiplication with same-center terms).
 *
 *  Typical usage
 *  -------------
 *  - Build analytic functions and evaluate them at given points (@ref evalf).
 *  - Construct @ref GaussExp<D> (expansions) by appending multiple @c GaussFunc.
 *  - Form products using @ref mult (returns @ref GaussPoly) or scale by a scalar.
 *  - Differentiate analytically with respect to a coordinate (@ref differentiate).
 */
template <int D> class GaussFunc : public Gaussian<D> {
public:
    /** @name Constructors
     *  @{
     */
    /** @brief Construct with isotropic exponent.
     *  @param beta   Isotropic exponent \f$\beta\f$ (used on all axes).
     *  @param alpha  Coefficient \f$\alpha\f$.
     *  @param pos    Center \f$\mathbf{R}\f$ (defaults to origin).
     *  @param pow    Powers \f$\mathbf{p}\f$ (defaults to all zeros).
     *
     *  This forwards to the @ref Gaussian<D> base constructor.
     */
    GaussFunc(double beta, double alpha, const Coord<D> &pos = {}, const std::array<int, D> &pow = {})
            : Gaussian<D>(beta, alpha, pos, pow) {}

    /** @brief Construct with anisotropic exponents (per-axis @p beta).
     *  @param beta   Array of exponents \f$(\beta_0,\dots,\beta_{D-1})\f$.
     *  @param alpha  Coefficient \f$\alpha\f$.
     *  @param pos    Center \f$\mathbf{R}\f$.
     *  @param pow    Powers \f$\mathbf{p}\f$.
     */
    GaussFunc(const std::array<double, D> &beta,
              double alpha,
              const Coord<D> &pos = {},
              const std::array<int, D> &pow = {})
            : Gaussian<D>(beta, alpha, pos, pow) {}

    /** @brief Copy constructor (shallow copy of POD members, as expected). */
    GaussFunc(const GaussFunc<D> &gf)
            : Gaussian<D>(gf) {}
    /** @brief Deleted assignment for safety (use copy-construct as needed). */
    GaussFunc<D> &operator=(const GaussFunc<D> &rhs) = delete;
    /** @brief Polymorphic copier (virtual constructor idiom). */
    Gaussian<D> *copy() const override;
    /** @} */

    /** @name Physics / analysis helpers
     *  @{
     */
    /** @brief Coulomb repulsion with another Gaussian (specialized for D=3).
     *  @details For D=3 and isotropic exponents, a closed form using Boys @f$F_0@f$ is used.
     *           Other D are not implemented and trigger a runtime abort.
     */
    double calcCoulombEnergy(const GaussFunc<D> &rhs) const;

    /** @brief \f$\|G\|_2^2 = \int |G|^2 \, d\mathbf{x}\f$ (separable product of 1D integrals). */
    double calcSquareNorm() const override;
    /** @} */

    /** @name Evaluation
     *  @{
     */
    /** @brief Full D-dimensional evaluation at coordinate @p r. */
    double evalf(const Coord<D> &r) const override;

    /** @brief 1D factor evaluation for axis @p dir (used in separable algorithms). */
    double evalf1D(double r, int dir) const override;
    /** @} */

    /** @name Transformations and algebra
     *  @{
     */
    /** @brief Wrap this single Gaussian as a length-1 Gaussian expansion. */
    GaussExp<D> asGaussExp() const override;

    /** @brief Analytic derivative w.r.t. @p dir, returns a @ref GaussPoly<D>. */
    GaussPoly<D> differentiate(int dir) const override;

    /** @brief In-place product with another Gaussian at the *same center*.
     *  @details Exponents and powers add; coefficients multiply.
     *           Fails fast if centers differ (cannot keep a pure GaussFunc).
     */
    void multInPlace(const GaussFunc<D> &rhs);
    /** @brief Alias for @ref multInPlace. */
    void operator*=(const GaussFunc<D> &rhs) { multInPlace(rhs); }

    /** @brief Product with another Gaussian (same or different center).
     *  @details Returns a @ref GaussPoly<D> (Gaussian times polynomial) obtained by
     *           completing the square and combining monomial factors. */
    GaussPoly<D> mult(const GaussFunc<D> &rhs);

    /** @brief Scalar multiplication (returns a scaled copy). */
    GaussFunc<D> mult(double c);

    /** @brief Operator overloads forwarding to the methods above. */
    GaussPoly<D> operator*(const GaussFunc<D> &rhs) { return this->mult(rhs); }
    GaussFunc<D> operator*(double c) { return this->mult(c); }
    /** @} */

    /** @name Power setters
     *  @{
     */
    /** @brief Set a single Cartesian power component @p power on axis @p d. */
    void setPow(int d, int power) override { this->power[d] = power; }
    /** @brief Set the full power vector \f$\mathbf{p}\f$. */
    void setPow(const std::array<int, D> &power) override { this->power = power; }
    /** @} */

private:
    /** @brief Pretty-printer used by stream insertion (see implementation). */
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp