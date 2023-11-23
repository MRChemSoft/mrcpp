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

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

namespace mrcpp {

/** @class JpowerIntegrals
 *
 * @brief A class needed for construction Schrodinger time evolution operator
 *
 * @details A two dimensional array consisting of integrals \f$ J_m \f$ as follows.
 * Our main operator has the following expansion
 * \f[
 *   \left[ \sigma_l^{\mathfrak n} \right]_{pj}
 *   (a)
 *   =
 *   \sum_{k = 0}^{\infty}
 *   C_{jp}^{2k}
 *   \widetilde J_{2k + j + p}(l, a)
 *   ,
 * \f]
 * where \f$ a = t \mathfrak N^2 = t 4^{\mathfrak n} \f$
 * and
 * \f[
 *     \widetilde J_m
 *     =
 *     \frac
 *     {
 *         I_m
 *         e^{ i \frac {\pi}4 (m - 1) }
 *     }
 *     {
 *         2 \pi ( m + 2 )!
 *     }
 *     =
 *     \frac
 *     {
 *         e^{ i \frac {\pi}4 (m - 1) }
 *     }
 *     {
 *         2 \pi ( m + 2 )!
 *     }
 *     \int_{\mathbb R}
 *     \exp
 *     \left(
 *         \rho l \exp \left( i \frac \pi 4 \right) - a \rho^2
 *     \right)
 *     \rho^m
 *     d \rho
 * \f]
 * satisfying the following relation
 * \f[
 *     \widetilde J_{m+1}
 *     =
 *     \frac
 *     {
 *         il
 *     }
 *     {
 *         2a (m + 3)
 *     }
 *     \widetilde J_m
 *     +
 *     \frac {im}{2a(m + 2)(m + 3)}
 *     \widetilde J_{m-1}
 *     =
 *     \frac
 *     {
 *         i
 *     }
 *     {
 *         2a (m + 3)
 *     }
 *     \left(
 *         l
 *         \widetilde J_m
 *         +
 *         \frac {m}{(m + 2)}
 *         \widetilde J_{m-1}
 *     \right)
 *     , \quad
 *     m = 0, 1, 2, \ldots,
 * \f]
 * with \f$ \widetilde J_{-1} = 0 \f$ and
 * \f[
 * \label{power_integral_0}
 *     \widetilde J_0
 *     =
 *     \frac{ e^{ -i \frac{\pi}4 } }{ 4 \sqrt{ \pi a } }
 *     \exp
 *     \left(
 *         \frac{il^2}{4a}
 *     \right)
 *     .
 * \f]
 *
 *  
 */
class JpowerIntegrals
{
public:
    /// @brief creates an array of power integrals
    /// @param a : parameter a
    /// @param scaling : scaling level
    /// @param M : maximum amount of integrals for each l
    /// @param threshold : lower limit for neglecting the integrals
    /// @details The array is orginised as a vector ordered as \f$l = 0, 1, 2, \ldots, 2^n - 1, 1 - 2^n, 2 - 2^n, \ldots, -2, -1 \f$.
    JpowerIntegrals(double a, int scaling, int M, double threshold = 1.0e-15);
    //JpowerIntegrals(const JpowerIntegrals& other);

    
    int scaling;  //it is probably not used
    std::vector<std::vector<std::complex<double>>> integrals;

    std::vector<std::complex<double>> & operator[](int index);
private:
    std::vector<std::complex<double>> calculate_J_power_integrals(int l, double a, int M, double threshold);
    void crop(std::vector<std::complex<double>> & J, double threshold);
};

} // namespace mrcpp
