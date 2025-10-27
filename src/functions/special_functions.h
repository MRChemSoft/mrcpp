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

#include <cmath>
#include <complex>

namespace mrcpp {

/**
 * # Free-particle Gaussian propagation (analytic form)
 *
 * @brief Analytic solution \f$\psi(x,t)\f$ of the 1D free-particle Schrödinger equation
 *        for a Gaussian initially centered at \f$x_0\f$ with width parameter \f$\sigma>0\f$.
 *
 * This declaration corresponds to the definition in `special_functions.cpp`. The solution used is
 * \f[
 *   \psi(x,t)
 *   =
 *   \sqrt{\frac{\sigma}{\,\sigma + 4\, i\, t\,}}\;
 *   \exp\!\left(-\,\frac{(x-x_0)^2}{\,\sigma + 4\, i\, t\,}\right),
 * \f]
 * which matches the conventional free propagator with units chosen such that \f$\hbar=1\f$
 * and mass \f$m=\tfrac12\f$ (hence the factor \f$4it\f$ in the denominator).
 *
 * @param x     Real-space coordinate \f$x \in \mathbb{R}\f$.
 * @param x0    Initial center \f$x_0\f$ of the Gaussian at \f$t=0\f$.
 * @param t     Time \f$t \in \mathbb{R}\f$.
 * @param sigma Positive width parameter \f$\sigma>0\f$ of the initial Gaussian.
 *
 * @return Complex value of \f$\psi(x,t)\f$ at the requested point.
 *
 * @note The complex square root in the prefactor is taken on the principal branch
 *       by `std::sqrt(std::complex<double>)`.
 * @note For \f$t=0\f$, the expression reduces to \f$\psi(x,0)=\exp\!\big(-\tfrac{(x-x_0)^2}{\sigma}\big)\f$.
 */
std::complex<double> free_particle_analytical_solution(double x, double x0, double t, double sigma);

/**
 * # Smooth compactly supported bump
 *
 * @brief A smooth (\f$C^\infty\f$) non-negative function supported on the open interval \f$(a,b)\f$.
 *
 * The function is defined by
 * \f[
 *   g_{a,b}(x) =
 *   \begin{cases}
 *     \exp\!\Big(-\dfrac{b-a}{(x-a)(b-x)}\Big), & a < x < b,\\[6pt]
 *     0, & \text{otherwise},
 *   \end{cases}
 * \f]
 * and vanishes to **all orders** at the endpoints \f$a\f$ and \f$b\f$.
 *
 * @param x Real-space coordinate \f$x \in \mathbb{R}\f$.
 * @param a Left endpoint of support (default `0`).
 * @param b Right endpoint of support (default `1`).
 *
 * @return The value \f$g_{a,b}(x)\f$.
 *
 * @note If \f$a \ge b\f$, the function is identically zero for all \f$x\f$.
 * @warning Near the endpoints, the denominator \f$(x-a)(b-x)\f$ becomes small;
 *          the exponent is large and negative so the result underflows smoothly to zero.
 */
double smooth_compact_function(double x, double a = 0, double b = 1);

} // namespace mrcpp