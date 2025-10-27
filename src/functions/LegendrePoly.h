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

#include "Polynomial.h"

namespace mrcpp {

/** @class LegendrePoly
 *  @brief Polynomial subclass representing a (possibly shifted/scaled) Legendre polynomial.
 *
 *  Purpose
 *  -------
 *  Encapsulates the Legendre polynomial \(P_k\) of degree @p k, constructed on the
 *  canonical interval \([-1,1]\) and then affinely mapped to an external coordinate
 *  via the Polynomial base class’ internal transform:
 *
 *    \f[
 *      q = N\,x + L,
 *    \f]
 *
 *  so that evaluations are effectively \(P_k(q(x))\).
 *
 *  Construction details
 *  --------------------
 *  - The raw coefficients of \(P_k(q)\) on \([-1,1]\) are computed using the
 *    standard three–term recurrence in @ref computeLegendrePolynomial.
 *  - After coefficients are set, the base class is instructed to translate by @p l
 *    and dilate by @p n, which stores the affine map \((N,L)\) used at evaluation time.
 *
 *  Notes
 *  -----
 *  - The actual caching of lower-order polynomials and the affine setup are handled
 *    in the corresponding .cpp file (see constructor and implementation comments).
 *  - Derivative helpers return values with respect to the *external* variable @p x,
 *    taking the internal affine map into account.
 */
class LegendrePoly final : public Polynomial {
public:
    /** @brief Construct degree-@p k Legendre polynomial with optional affine transform.
     *
     *  @param k Degree (order) of the Legendre polynomial \(P_k\).
     *  @param n Dilation factor (applied after translation). Conceptually produces \(P_k(Nx+L)\) with \(N=n\).
     *  @param l Translation (applied before dilation). Conceptually produces \(P_k(Nx+L)\) with \(L=l\).
     *
     *  Semantics
     *  ---------
     *  - First builds \(P_k\) on \([-1,1]\) in the internal variable \(q\).
     *  - Records canonical bounds \([-1,1]\) for error checking.
     *  - Applies the affine map encoded by @p n and @p l through the base class.
     */
    LegendrePoly(int k, double n = 1.0, double l = 0.0);

    /** @brief Evaluate \(P_k(x)\) and its first derivative w.r.t. the external variable.
     *
     *  @param x External evaluation point.
     *  @return \f$[\,P_k(x),\,\tfrac{d}{dx}P_k(x)\,]\f$ as an Eigen::Vector2d.
     *
     *  Details
     *  -------
     *  - Internally maps @p x to the polynomial’s canonical coordinate \(q = N x + L\).
     *  - Uses a recurrence that simultaneously advances value and derivative.
     *  - Performs a bounds check consistent with the base-class domain bookkeeping.
     */
    Eigen::Vector2d firstDerivative(double x) const;

    /** @brief Evaluate value, first and second derivatives (declared interface).
     *
     *  @param x External evaluation point.
     *  @return \f$[\,P_k(x),\,P'_k(x),\,P''_k(x)\,]\f$ as an Eigen::Vector3d.
     *
     *  @note The current implementation in the .cpp intentionally aborts
     *        (NOT_IMPLEMENTED) to document that second-derivative support
     *        is not provided yet.
     */
    Eigen::Vector3d secondDerivative(double x) const;

private:
    /** @brief Fill coefficient vector with the canonical \([-1,1]\) Legendre polynomial \(P_k\).
     *
     *  @param k Degree (order).
     *
     *  Implementation sketch
     *  ---------------------
     *  - Base cases:
     *      - \(P_0(q) = 1\)
     *      - \(P_1(q) = q\)
     *  - Recurrence for \(k \ge 2\):
     *      \f[
     *        P_k(q) = \frac{(2k-1)\,q\,P_{k-1}(q) - (k-1)\,P_{k-2}(q)}{k}.
     *      \f]
     *  - Operates directly in coefficient space (ascending powers of \(q\)).
     *  - Lower orders are retrieved from an ObjectCache to avoid recomputation (in .cpp).
     */
    void computeLegendrePolynomial(int k);
};

} // namespace mrcpp