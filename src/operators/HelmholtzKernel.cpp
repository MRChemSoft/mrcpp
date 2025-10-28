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
 * @file HelmholtzKernel.cpp
 * @brief Gaussian expansion of the 3D screened Coulomb / Helmholtz kernel.
 *
 * @details
 * This file implements a separable Gaussian approximation to the radial 3D
 * Helmholtz kernel on a finite interval \f$[r_\text{min}, r_\text{max}]\f$.
 * The expansion has the form
 * \f[
 *   K_\mu(r) \;\approx\; \sum_{m=1}^{M} \beta_m\, e^{-\alpha_m r^2},
 * \f]
 * where \f$M\f$ (the separation rank) and the parameters \f$\{\alpha_m,\beta_m\}\f$
 * are chosen by truncating and sampling an integral representation with a uniform
 * trapezoidal rule in a logarithmic variable \f$s\f$. The resulting coefficients
 * depend on:
 *  - the screening parameter \f$\mu > 0\f$,
 *  - a target relative accuracy \f$\varepsilon\f$,
 *  - a radial domain \f$[r_\text{min}, r_\text{max}]\f$.
 *
 * Internally, the interval is rescaled to \f$[r_\text{min}/r_\text{max}, 1]\f$
 * to keep the step-size heuristics well-conditioned; the generated Gaussian
 * parameters are then rescaled back to the original units.
 */

#include "HelmholtzKernel.h"

#include <cmath>

#include "functions/GaussFunc.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @class HelmholtzKernel
 * @brief Builds a 1D Gaussian expansion that approximates the 3D Helmholtz kernel.
 *
 * @details
 * The constructor discretizes an auxiliary integral over a log-scaled variable
 * \f$s\fin[s_1,s_2]\f$ using a uniform step \f$h\f$ derived from the requested
 * tolerance \f$\varepsilon\f$. For each quadrature node it produces a single
 * Gaussian term with exponent \f$\alpha_m\f$ and weight \f$\beta_m\f$. Endpoints
 * receive the trapezoidal half-weights.
 *
 * Rescaling:
 * - Define \f$r_0 = r_\text{min}/r_\text{max}\f$ and \f$r_1 = r_\text{max}\f$.
 * - Work on \f$[r_0,1]\f$, then map back by multiplying
 *   \f$\alpha \leftarrow \alpha / r_1^2\f$ and \f$\beta \leftarrow \beta / r_1\f$.
 *
 * Rank control:
 * - The number of exponentials is \f$M = \lceil (s_2 - s_1)/h \rceil + 1\f$.
 * - If \f$M > \texttt{MaxSepRank}\f$ the constructor aborts, signaling that the
 *   requested accuracy on the given domain would require too large a rank.
 *
 * @param mu       Screening parameter \f$\mu > 0\f$.
 * @param epsilon  Target relative accuracy \f$\varepsilon \in (0,1)\f$.
 * @param r_min    Minimal radius of the approximation interval (strictly positive).
 * @param r_max    Maximal radius of the approximation interval (\f$r_\text{max} > r_\text{min}\f$).
 *
 * @note
 * This routine assumes the standard MRCPP constants \c pi and \c root_pi are available
 * in the \c mrcpp namespace and that \c MaxSepRank bounds the admissible separation rank.
 */
HelmholtzKernel::HelmholtzKernel(double mu, double epsilon, double r_min, double r_max)
        : GaussExp<1>() {
    // Rescale the interval to [r0, 1] and precompute scaled mu
    const double r0 = r_min / r_max;
    const double r1 = r_max;
    const double mu_tilde = mu * r1;

    // Truncation window [s1, s2] giving ~epsilon relative error
    // The heuristic t = max(-2.5 ln eps, 5) balances tails for practical eps
    const long double t = std::max((-2.5L * std::log(epsilon)), 5.0L);
    const double s1 = -std::log(4.0L * t / (mu_tilde * mu_tilde)) / 2.0L;
    const double s2 =  std::log(t / (r0 * r0)) / 2.0L;

    // Trapezoidal step size h from an empirical fit versus log10(epsilon)
    const double h = 1.0 / (0.20L - 0.47L * std::log10(epsilon));
    const int n_exp = static_cast<int>(std::ceil((s2 - s1) / h) + 1.0);

    if (n_exp > MaxSepRank) MSG_ABORT("Maximum separation rank exceeded.");

    // Uniform trapezoidal quadrature in s; endpoints get half-weight.
    for (int i = 0; i < n_exp; ++i) {
        const double s = s1 + h * i;

        // Intermediate quantities (written explicitly for clarity)
        // temp  = -2 s
        // temp2 = - (mu_tilde^2) e^{-2 s} / 4 + s
        // beta  ~ h * 2/sqrt(pi) * exp(temp2)
        const double temp  = -2.0 * s;
        const double temp2 = - (mu_tilde * mu_tilde) * std::exp(temp) / 4.0 + s;

        double beta  = h * (2.0 / root_pi) * std::exp(temp2);
        double alpha = std::exp(2.0L * s);

        // Rescale back to the original radial units
        alpha *= 1.0 / (r1 * r1);
        beta  *= 1.0 / r1;

        // Trapezoidal half-weights at the endpoints
        if (i == 0 || i == (n_exp - 1)) beta *= 0.5;

        // Append the 1D Gaussian term exp(-alpha r^2) with prefactor beta
        GaussFunc<1> gFunc(alpha, beta);
        this->append(gFunc);
    }
}

} // namespace mrcpp