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


/** @brief Free-particle time evolution on real line.
 *
 * @param[in] x: space coordinate in \f$ \mathbb R \f$.
 * @param[in] x0: \f$ x_0 \f$ center of gaussian function at zero time moment.
 * @param[in] t: time moment.
 * @param[in] sigma: \f$ \sigma \f$ width of the initial gaussian wave.
 *
 * @details Analytical solution of a one dimensional free-particle
 * movement
 * \f[
 *      \psi(x, t)
 *      =
 *      \sqrt{
 *          \frac{ \sigma }{ 4it + \sigma }
 *      }
 *      e^{ - \frac { (x - x_0)^2 }{ 4it + \sigma } }
 * \f]
 * where \f$ t, \sigma > 0 \f$.
 * 
 * @returns The complex-valued wave function
 * \f$ \psi(x, t) \f$
 * at the specified space coordinate and time.
 * 
 * 
 */
std::complex<double> free_particle_analytical_solution(double x, double x0, double t, double sigma)
{
    std::complex<double> i(0.0, 1.0);  // Imaginary unit
    auto denominator = 4 * t * i + sigma;
    std::complex<double> sqrt_denom = std::sqrt(denominator);
    std::complex<double> exponent = -((x - x0) * (x - x0)) / denominator;

    return std::sqrt(sigma) / sqrt_denom * std::exp(exponent);
}



/** @brief A smooth compactly supported non-negative function.
 *
 * @param[in] x: space coordinate in \f$ \mathbb R \f$.
 * @param[in] a: the left support boundary.
 * @param[in] b: the right support boundary.
 *
 * @details Smooth function on the real line \f$ \mathbb R \f$
 * defined by the formula
 * \f[
 *      g_{a,b} (x) = \exp \left( - \frac{b - a}{(x - a)(b - x)} \right)
 *      , \quad
 *      a < x < b
 * \f]
 * and \f$ g_{a,b} (x) = 0 \f$ elsewhere.
 * 
 * @returns The non-negative value
 * \f$ g_{a,b} (x) \f$
 * at the specified space coordinate \f$ x \in \mathbb R \f$.
 * 
 * 
 */
double smooth_compact_function(double x, double a, double b) {
    double res = 0;
    if (a < x && x < b) {
        res = exp((a - b) / (x - a) / (b - x));
    }
    return res;
}

} // namespace mrcpp