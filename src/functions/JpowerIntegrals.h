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
 * @brief A class needed for construction of the Schrödinger time-evolution operator
 */
class JpowerIntegrals {
public:
    /// Create the table of power integrals
    /// The array is organised as l = 0,1,...,2^n-1, 1-2^n, 2-2^n, ..., -2, -1
    JpowerIntegrals(double a, int scaling, int M, double threshold = 1.0e-15);

    int scaling;  // may be unused
    std::vector<std::vector<std::complex<double>>> integrals;

    std::vector<std::complex<double>> &operator[](int index);

private:
    std::vector<std::complex<double>> calculate_J_power_integrals(int l, double a, int M, double threshold);
    void crop(std::vector<std::complex<double>> &J, double threshold);
};

/** @class DerivativePowerIntegrals
 *
 * @brief Power-integral table used by the smooth-derivative operator.
 *
 * Implementation uses an FFT-based construction in the .cpp (FFTW is included
 * there under MRCPP_HAVE_FFTW). Header intentionally does not include fftw3.h.
 */
class DerivativePowerIntegrals {
public:
    /// Build the table for a given scaling level and cutoff
    DerivativePowerIntegrals(double cut_off, int scaling, int M, double threshold = 1.0e-15);

    int scaling;
    // same l-indexing convention as JpowerIntegrals; entries are real-valued
    std::vector<std::vector<double>> integrals;

    // Access by l in [-2^n+1, ..., 2^n-1] (negative l mapped to the tail)
    std::vector<double> &operator[](int index);

private:
    std::vector<std::vector<double>> calculate_J_power_integrals(double cut_off, int M, double threshold);
};

} // namespace mrcpp