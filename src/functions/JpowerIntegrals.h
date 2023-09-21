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

class JpowerIntegrals
{
public:
    /// @brief creates an array of power integrals
    /// @param a : parameter a
    /// @param N : 2^n
    /// @param M : maximum amount of integrals for each l
    /// @param treshold : lower limit for neglecting the integrals
    JpowerIntegrals(double a, int N, int M, double treshold = 1.0e-15);

    ~JpowerIntegrals();

    std::vector<std::vector<std::complex<double>>> integrals;

    std::vector<std::complex<double>> & operator[](int index);
    std::vector<std::complex<double>> calculate_J_power_inetgarls(int l, double a, int M, double treshold);
};

} // namespace mrcpp
