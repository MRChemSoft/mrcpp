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
 * @file PoissonKernel.cpp
 * @brief Builds a Gaussian expansion approximation of the 3D Poisson kernel.
 *
 * @details
 * This implementation constructs a separated approximation to the radial
 * Poisson kernel
 * \f[
 *   \frac{1}{\lvert \mathbf r \rvert}
 * \f]
 * on a finite annulus \f$ r \in [r_{\min},\, r_{\max}] \f$ by means of a
 * finite sum of Gaussians
 * \f[
 *   \frac{1}{r} \;\approx\; \sum_{m=1}^{M} \beta_m \, e^{-\alpha_m r^2},
 * \f]
 * where the coefficients \f$ \{\alpha_m,\beta_m\} \f$ are obtained by
 * truncating and discretizing (via the trapezoidal rule) a suitable integral
 * representation of \f$ 1/r \f$ in logarithmic variables. The truncation
 * bounds \f$[s_1, s_2]\f$ and the step \f$h\f$ are chosen to meet a requested
 * relative accuracy \c epsilon on the normalized interval \f$[r_{\min}/r_{\max},\,1]\f$,
 * after which the expansion is rescaled back to \f$[r_{\min},\,r_{\max}]\f$.
 *
 * ### Inputs
 * - \c epsilon: Target relative error for the expansion (heuristic, affects
 *   the truncation window and step size).
 * - \c r_min, \c r_max: Inner/outer radii that define the interval of validity.
 *
 * ### Algorithm sketch
 * 1. Normalize the domain to \f$[r_0, 1]\f$ with \f$r_0 = r_{\min}/r_{\max}\f$ and set
 *    \f$r_1 = r_{\max}\f$ for subsequent rescaling.
 * 2. Determine auxiliary parameters \f$t_1, t_2\f$ such that the tails of the
 *    integral representation are below \c epsilon.
 * 3. Convert tails to truncation limits \f$s_1, s_2\f$ in logarithmic coordinates.
 * 4. Choose trapezoidal step size \f$h\f$ as a function of \c epsilon and compute
 *    the number of terms \f$M\f$.
 * 5. Form nodes \f$s_i = s_1 + i h\f$ and corresponding Gaussian parameters
 *    \f$\alpha_i, \beta_i\f$ (with endpoint halving for the trapezoid rule).
 * 6. Rescale \f$\alpha_i, \beta_i\f$ from the normalized interval back to
 *    \f$[r_{\min}, r_{\max}]\f$ and append each term to the @ref GaussExp.
 *
 * The resulting expansion length is capped by \c MaxSepRank; exceeding this
 * limit aborts construction.
 */

#include "PoissonKernel.h"

#include <cmath>

#include "functions/GaussFunc.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @brief Construct a Gaussian expansion of the 3D Poisson kernel on \f$[r_{\min}, r_{\max}]\f$.
 *
 * @param epsilon Target relative accuracy for the expansion (heuristic).
 * @param r_min   Minimum radius of the interval of validity (\f$>0\f$).
 * @param r_max   Maximum radius of the interval of validity (\f$> r_{\min}\f$).
 *
 * @details
 * The method chooses truncation limits \f$s_1, s_2\f$ and a step size \f$h\f$
 * for a trapezoidal discretization so that the contribution of neglected tails
 * is below \c epsilon in the normalized variable. Each quadrature node yields
 * one Gaussian term. Endpoint weights are halved, as per the trapezoidal rule.
 *
 * The final expansion is rescaled to the physical interval by the mappings
 * \f$ \alpha \leftarrow \alpha / r_{\max}^2 \f$ and \f$ \beta \leftarrow \beta / r_{\max} \f$,
 * ensuring that the approximation targets the original (unscaled) radius.
 *
 * @note If the number of terms exceeds @c MaxSepRank, construction aborts.
 */
PoissonKernel::PoissonKernel(double epsilon, double r_min, double r_max)
        : GaussExp<1>() {
    // Constructed on [rMin/rMax, 1.0], then rescaled to [rMin, rMax]
    double r0 = r_min / r_max;
    double r1 = r_max;

    // Choose t1, t2 so that tail contributions are below epsilon
    double t1 = 1.0L;
    while ((2.0 * t1 * std::exp(-t1)) > epsilon) t1 *= 1.1L;

    double t2 = 1.0L;
    while ((std::sqrt(t2) * std::exp(-t2) / r0) > epsilon) t2 *= 1.1L;

    // Truncation window [s1, s2] ensuring relative error ~ epsilon
    double s1 = -std::log(2.0 * t1);
    double s2 = std::log(t2 / (r0 * r0)) / 2.0;

    // Trapezoidal step size h determined from epsilon (empirical fit)
    double h = 1.0 / (0.2L - 0.47L * std::log10(epsilon));
    int n_exp = static_cast<int>(std::ceil((s2 - s1) / h) + 1);
    if (n_exp > MaxSepRank) MSG_ABORT("Maximum separation rank exceeded.");

    for (int i = 0; i < n_exp; i++) {
        double arg = s1 + h * i;
        double sinharg = std::sinh(arg);
        double cosharg = std::cosh(arg);
        double onepexp = 1.0 + std::exp(-sinharg);

        // Parameters before rescaling back to [r_min, r_max]
        double expo = 4.0L * (sinharg + std::log(onepexp)) * (sinharg + std::log(onepexp));
        double coef = h * (4.0L / root_pi) * cosharg / onepexp;

        // Rescale to physical interval
        expo *= 1.0 / (r1 * r1);
        coef *= 1.0 / r1;

        // Trapezoidal rule endpoint correction
        if (i == 0 || i == (n_exp - 1)) coef *= 1.0 / 2.0;

        GaussFunc<1> gFunc(expo, coef);
        this->append(gFunc);
    }
}

} // namespace mrcpp
