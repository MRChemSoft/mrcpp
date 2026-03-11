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
 * @file PoissonKernel.h
 * @brief Declaration of a Gaussian-expansion approximation to the 3D Poisson kernel.
 */

#pragma once

#include "functions/GaussExp.h"

namespace mrcpp {

/**
 * @class PoissonKernel
 * @brief Gaussian expansion of the radial Poisson kernel \f$ 1/r \f$ on a bounded interval.
 *
 * @details
 * Builds a separated, finite Gaussian expansion that approximates the 3D Poisson kernel
 * \f[
 *   \frac{1}{\lvert \mathbf r \rvert} \;\approx\; \sum_{m=1}^{M} \beta_m \, e^{-\alpha_m r^2},
 * \f]
 * valid for radii \f$ r \in [r_{\min},\, r_{\max}] \f$. The coefficients
 * \f$ \{\alpha_m,\beta_m\}_{m=1}^M \f$ are produced by truncating and discretizing
 * (via a trapezoidal rule in logarithmic variables) a continuous representation of
 * \f$ 1/r \f$, with the truncation window and step size chosen to meet a target
 * relative tolerance \p epsilon on the *normalized* interval
 * \f$ [r_{\min}/r_{\max},\,1] \f$ and then rescaled back to \f$ [r_{\min}, r_{\max}] \f$.
 *
 * The resulting object is a one-dimensional @ref GaussExp "GaussExp<1>" whose entries
 * can be used by separable convolution operators to assemble higher-dimensional
 * kernels and operators.
 *
 * @note
 * - Requires \f$ r_{\min} > 0 \f$ and \f$ r_{\max} > r_{\min} \f$.
 * - The number of Gaussian terms \f$ M \f$ is bounded internally (see `MaxSepRank`);
 *   exceeding this bound will abort construction in the implementation.
 *
 * @see GaussExp, GaussFunc
 */
class PoissonKernel final : public GaussExp<1> {
public:
    /**
     * @brief Construct a Gaussian expansion of \f$ 1/r \f$ on \f$ [r_{\min}, r_{\max}] \f$.
     *
     * @param epsilon Target relative accuracy (heuristic; smaller ⇒ more terms).
     * @param r_min   Lower radius of validity, must satisfy \f$ r_{\min} > 0 \f$.
     * @param r_max   Upper radius of validity, must satisfy \f$ r_{\max} > r_{\min} \f$.
     *
     * @details
     * Populates this @ref GaussExp with terms \f$ (\alpha_m,\beta_m) \f$ so that
     * \f$ \sum_m \beta_m e^{-\alpha_m r^2} \approx 1/r \f$ over the requested interval.
     * Coefficients are ordered according to the underlying quadrature and include
     * standard endpoint weighting for the trapezoidal rule.
     */
    PoissonKernel(double epsilon, double r_min, double r_max);
};

} // namespace mrcpp