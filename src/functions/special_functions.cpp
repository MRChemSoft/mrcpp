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

std::complex<double> free_particle_analytical_solution(double x, double x0, double t, double sigma)
{
    std::complex<double> i(0.0, 1.0);
    std::complex<double> denom = sigma + 4.0 * t * i;
    std::complex<double> exponent = -((x - x0) * (x - x0)) / denom;

    return std::sqrt(sigma) / std::sqrt(denom) * std::exp(exponent);
}

double smooth_compact_function(double x, double a, double b) {
    if (a < x && x < b) {
        return std::exp((a - b) / ((x - a) * (b - x)));
    }
    return 0.0;
}

} // namespace mrcpp