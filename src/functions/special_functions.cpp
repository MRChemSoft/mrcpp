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

#include "special_functions.h"

namespace mrcpp {

/**
 * @brief Analytic solution of the free-particle Schrödinger equation on ℝ at time @p t.
 *
 * This implements the standard Gaussian wave packet propagation (free particle, \f$\hbar=1\f$, mass \f$m=\tfrac12\f$
 * so that the free propagator denominator becomes \f$4it+\sigma\f$ as used below). Given an initial
 * Gaussian of width parameter \f$\sigma>0\f$ centered at \f$x_0\f$ at time \f$t=0\f$,
 * the wave function at time \f$t\f$ is
 *
 * \f[
 *   \psi(x,t)
 *   =
 *   \sqrt{\frac{\sigma}{\,\sigma + 4\, i\, t\,}}
 *   \exp\!\left(
 *     -\,\frac{(x - x_0)^2}{\,\sigma + 4\, i\, t\,}
 *   \right),
 * \f]
 *
 * which disperses in time and acquires a complex phase.
 *
 * #### Parameters
 * - @param x     Real-space coordinate \f$x \in \mathbb{R}\f$.
 * - @param x0    Initial center \f$x_0\f$ of the Gaussian at \f$t=0\f$.
 * - @param t     Time \f$t \in \mathbb{R}\f$ (can be positive or negative).
 * - @param sigma Width parameter \f$\sigma>0\f$ of the initial Gaussian.
 *
 * #### Returns
 * The complex-valued wave function \f$\psi(x,t)\f$ at the requested space-time point.
 *
 * #### Notes
 * - For @p t = 0, this reduces to \f$\psi(x,0)=\exp\!\big(-\tfrac{(x-x_0)^2}{\sigma}\big)\f$.
 * - The branch of the complex square root is the principal branch via `std::sqrt(std::complex)`.
 * - Numerical behavior near large \f$|t|\f$: the modulus decays like \f$|\sigma/(\sigma+4it)|^{1/2}\f$,
 *   while the phase is dominated by the complex denominator; standard `std::complex` arithmetic handles this.
 * - This function assumes consistent physical units so that the closed form above applies directly.
 */
std::complex<double> free_particle_analytical_solution(double x, double x0, double t, double sigma)
{
    std::complex<double> i(0.0, 1.0);                  // imaginary unit i
    std::complex<double> denom = sigma + 4.0 * t * i;  // σ + 4 i t
    std::complex<double> exponent = -((x - x0) * (x - x0)) / denom;

    return std::sqrt(sigma) / std::sqrt(denom) * std::exp(exponent);
}

/**
 * @brief Smooth, compactly supported "bump" function on the interval \f$(a,b)\f$.
 *
 * Defines a non-negative \f$C^\infty\f$ function
 * \f[
 *   g_{a,b}(x) =
 *   \begin{cases}
 *     \exp\!\Big( -\,\dfrac{b-a}{(x-a)(b-x)} \Big), & a < x < b,\\[6pt]
 *     0, & \text{otherwise},
 *   \end{cases}
 * \f]
 * which vanishes to **all orders** at the endpoints \f$a\f$ and \f$b\f$.
 *
 * #### Parameters
 * - @param x  Real-space coordinate \f$x \in \mathbb{R}\f$.
 * - @param a  Left endpoint (must satisfy \f$a<b\f$ for a non-trivial function).
 * - @param b  Right endpoint (\f$b>a\f$).
 *
 * #### Returns
 * - \f$g_{a,b}(x)\f$ if \f$a < x < b\f$, and `0.0` otherwise.
 *
 * #### Numerical remarks
 * - Near the endpoints, \f$(x-a)(b-x)\to 0^+\f$ and the exponent \f$-\frac{b-a}{(x-a)(b-x)}\f$ becomes large
 *   and negative, so the value safely underflows toward 0; this is expected and preserves smooth compact support.
 * - If `a >= b`, the definition yields the zero function for all `x`.
 */
double smooth_compact_function(double x, double a, double b) {
    if (a < x && x < b) {
        // Equivalent to: exp( - (b-a) / ((x-a)(b-x)) )
        return std::exp((a - b) / ((x - a) * (b - x)));
    }
    return 0.0;
}

} // namespace mrcpp