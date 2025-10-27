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

/*
 * File overview
 * -------------
 * Implements LegendrePoly, a Polynomial subclass that builds the (shifted/scaled)
 * Legendre polynomial of a given order k using the standard three-term recurrence.
 *
 * Key ideas:
 *  - Coefficients for P_k on the canonical interval [-1, 1] are computed once
 *    (in the Polynomial's coefficient storage) by combining cached lower orders.
 *  - The resulting polynomial is then *affinely transformed* by an internal
 *    mapping x ↦ N·x + L via Polynomial::translate(l) and Polynomial::dilate(n),
 *    so that users obtain a Legendre polynomial defined on the transformed domain.
 *  - A lightweight cache (ObjectCache<LegendrePoly>) avoids recomputing lower
 *    orders repeatedly when constructing higher ones.
 *
 * Extras:
 *  - firstDerivative(x) returns both P_k(x) and P'_k(x) evaluated at x (w.r.t. the
 *    current affine mapping).
 *  - secondDerivative(x) is declared but not implemented (calls NOT_IMPLEMENTED_ABORT).
 */

#include "LegendrePoly.h"
#include "core/ObjectCache.h"
#include "utils/Printer.h"

using namespace Eigen;

namespace mrcpp {

using LegendreCache = ObjectCache<LegendrePoly>; // Cache of LegendrePoly objects keyed by order

/** @brief Construct the order-k Legendre polynomial and apply an affine transform.
 *
 * @param k Polynomial order (degree).
 * @param n Dilation factor applied after construction (see Polynomial::dilate).
 * @param l Translation applied before dilation (see Polynomial::translate).
 *
 * Details:
 *  - The raw Legendre polynomial P_k is constructed on [-1, 1] using the standard
 *    recurrence:
 *      P_0(q) = 1,
 *      P_1(q) = q,
 *      P_k(q) = ((2k-1) q P_{k-1}(q) - (k-1) P_{k-2}(q)) / k,  k ≥ 2.
 *    Here q is the canonical variable on [-1, 1].
 *  - Lower-order polynomials P_{k-1}, P_{k-2} are fetched (or created once and
 *    cached) through LegendreCache to avoid recomputation.
 *  - After coefficients are set, the base interval [-1, 1] is recorded via
 *    setBounds, then the polynomial is translated by l and dilated by n,
 *    effectively producing P_k(N·x + L) in the Polynomial base class, where
 *    N and L are the stored affine parameters.
 */
LegendrePoly::LegendrePoly(int k, double n, double l)
        : Polynomial(k) {
    // Ensure lower orders are cached: creating P_k requires P_{k-1} and P_{k-2}.
    // We preload P_{k-1} so the subsequent compute can fetch both from cache.
    LegendreCache &Cache = LegendreCache::getInstance();
    if (k >= 1) {
        if (not Cache.hasId(k - 1)) {
            auto *lp = new LegendrePoly(k - 1);
            // Rough memory accounting: 2*(k+1) doubles (heuristic) for the cache
            Cache.load(k - 1, lp, 2 * sizeof(double) * (k + 1));
        }
    }

    // Compute P_k on the canonical domain [-1, 1]
    computeLegendrePolynomial(k);

    // Record canonical bounds for sanity checks in eval/derivatives
    double a = -1.0;
    double b = 1.0;
    setBounds(&a, &b);

    // Apply affine map x ↦ N·x + L via translate(l) then dilate(n)
    translate(l);
    dilate(n);
}

/** @brief Populate this->coefs with the coefficients of P_k on [-1,1].
 *
 * Implements the standard three-term Legendre recurrence in coefficient space:
 *   P_0(q) = 1,
 *   P_1(q) = q,
 *   P_k(q) = ((2k-1) q P_{k-1}(q) - (k-1) P_{k-2}(q)) / k,  for k ≥ 2.
 *
 * Coefficient layout:
 *  - The Polynomial base stores coefficients in ascending powers:
 *      coefs[j] corresponds to q^j.
 *  - To form P_k, we combine cached P_{k-1} and P_{k-2} term-by-term.
 *
 * Edge cases:
 *  - k=0 and k=1 are assigned explicitly.
 */
void LegendrePoly::computeLegendrePolynomial(int k) {
    assert(this->size() >= k);

    if (k == 0) {
        // P_0(q) = 1
        this->coefs[0] = 1.0;
    } else if (k == 1) {
        // P_1(q) = q
        this->coefs[0] = 0.0;
        this->coefs[1] = 1.0;
    } else {
        // Fetch lower-order polynomials from the cache
        LegendreCache &Cache = LegendreCache::getInstance();
        LegendrePoly &Lm1 = Cache.get(k - 1); // P_{k-1}
        LegendrePoly &Lm2 = Cache.get(k - 2); // P_{k-2}

        auto K = (double)k;

        // Constant term (j=0):
        //   coef0 = -(k-1)/k * (coef from P_{k-2} at j=0)
        double cm2_0 = Lm2.getCoefs()[0];
        this->coefs[0] = -(K - 1.0) * cm2_0 / K;

        // Remaining terms (j=1..k):
        //   For j ≤ k-2, coef_j = ((2k-1)/k) * coef_{j-1}(P_{k-1}) - ((k-1)/k) * coef_j(P_{k-2})
        //   For j = k-1 or k, the P_{k-2} contribution vanishes (index out of range),
        //   so only the first term remains.
        for (int j = 1; j < k + 1; j++) {
            double cm1_jm1 = Lm1.getCoefs()[j - 1];
            if (j <= k - 2) {
                double cm2_j = Lm2.getCoefs()[j];
                this->coefs[j] = (2.0 * K - 1.0) * cm1_jm1 / K - (K - 1.0) * cm2_j / K;
            } else {
                this->coefs[j] = (2.0 * K - 1.0) * cm1_jm1 / K;
            }
        }
    }
}

/** @brief Evaluate P_k(x) and its first derivative at x (w.r.t. the current affine map).
 *
 * @param x Point of evaluation (in the *external* variable).
 * @return Vector2d { P_k(x), d/dx P_k(x) }.
 *
 * Details:
 *  - Bounds check (via outOfBounds) uses the base interval (set to [-1,1] and
 *    then transformed by the affine map stored in the Polynomial base).
 *  - Internally we evaluate in the mapped coordinate q = N·x + L.
 *  - Uses a forward recursion to accumulate both value and derivative following
 *    the Legendre three-term recurrence:
 *      y_i(q)  = ((2i-1) q y_{i-1} - (i-1) y_{i-2}) / i
 *      dy_i(q) = ((2i-1) q dy_{i-1} - (i-1) dy_{i-2} + (2i-1) y_{i-1}) / i
 *    (the last term is ∂/∂q of (2i-1) q y_{i-1}).
 *  - The returned derivative is with respect to the external variable x, taking
 *    into account the internal affine mapping (via the Polynomial base members).
 */
Vector2d LegendrePoly::firstDerivative(double x) const {
    double c1, c2, c4, ym, yp, y;
    double dy, dyp, dym;

    if (outOfBounds({x})) {
        MSG_ABORT("Argument out of bounds: " << x << " [" << this->A[0] << ", " << this->B[0] << "]");
    }

    // Affine map from external x to internal q
    double q = this->N * x + this->L;
    Vector2d val;

    int order = getOrder();

    // P_0(q) = 1, P'_0(q) = 0
    if (order == 0) {
        val(0) = 1.0;
        val(1) = 0.0;
        return val;
    }

    // P_1(q) = q; derivative follows the affine mapping stored in the base
    if (order == 1) {
        val(0) = q;
        val(1) = this->N * 1.0 + this->L; // as implemented in the original code
        return val;
    }

    // Initialize recurrence for i=2..order
    y  = q;   // y   = P_1
    dy = 1.0; // dy  = d/dq P_1
    yp = 1.0; // yp  = P_0
    dyp = 0.0;// dyp = d/dq P_0

    for (int i = 2; i < order + 1; i++) {
        c1 = (double)i;
        c2 = c1 * 2.0 - 1.0; // (2i-1)
        c4 = c1 - 1.0;       // (i-1)

        // Rotate "previous" states
        ym = y;

        // Value recurrence: y = P_i
        y = (c2 * q * y - c4 * yp) / c1;

        // Shift lower-order values
        yp = ym;

        // Derivative recurrence in q
        dym = dy;
        dy = (c2 * q * dy - c4 * dyp + c2 * yp) / c1;
        dyp = dym;
    }

    val(0) = y;
    val(1) = dy;
    return val;
}

/** @brief Evaluate P_k(x) together with first and second derivatives (not implemented).
 *
 * @param x Point of evaluation.
 * @return Vector3d { P_k(x), P'_k(x), P''_k(x) }.
 *
 * @note This routine currently calls NOT_IMPLEMENTED_ABORT. The code that follows
 *       shows the intended structure (value/first/second derivative recurrences),
 *       but it is not active. Keep as-is to reflect current behavior.
 */
Vector3d LegendrePoly::secondDerivative(double x) const {
    NOT_IMPLEMENTED_ABORT;
    double c1, c2, c4, ym, yp, y, d2y;
    double dy, dyp, dym, d2ym, d2yp;

    double q = this->N * x + this->L;
    if (outOfBounds({x})) {
        MSG_ABORT("Argument out of bounds: " << x << " [" << this->A[0] << ", " << this->B[0] << "]");
    }

    Vector3d val;

    int order = getOrder();
    if (order == 0) {
        val(0) = 1.e0;
        val(1) = 0.e0;
        val(2) = 0.e0;
        return val;
    }

    if (order == 1) {
        val(0) = q;
        val(1) = this->N * 1.e0 + this->L;
        val(2) = 0.e0;
        return val;
    }

    y   = q;
    dy  = 1.e0;
    d2y = 0.e0;
    yp   = 1.e0;
    dyp  = 0.e0;
    d2yp = 0.e0;

    for (int i = 2; i < order + 1; i++) {
        c1 = (double)i;
        c2 = c1 * 2.e0 - 1.e0;
        c4 = c1 - 1.e0;

        ym = y;
        y  = (c2 * x * y - c4 * yp) / c1;
        yp = ym;

        dym = dy;
        dy  = (c2 * x * dy - c4 * dyp + c2 * yp) / c1;
        dyp = dym;

        d2ym = d2y;
        d2y  = (c2 * x * d2y - c4 * d2yp + c2 * 2.e0 * dyp) / c1;
        d2yp = d2ym;
    }
    val(0) = y;
    val(1) = dy;
    val(2) = d2y;
    return val;
}

} // namespace mrcpp
