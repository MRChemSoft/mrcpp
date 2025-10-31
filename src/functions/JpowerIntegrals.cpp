/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3 of the License, or
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

#include "JpowerIntegrals.h"
#include <algorithm>

namespace mrcpp {

JpowerIntegrals::JpowerIntegrals(double a, int scaling, int M, double threshold) {
    this->scaling = scaling;
    int N = 1 << scaling;
    for (int l = 0; l < N; l++) integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
    for (int l = 1 - N; l < 0; l++) integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
}

std::vector<std::complex<double>> &JpowerIntegrals::operator[](int index) {
    if (index < 0) index += integrals.size();
    return integrals[index];
}

std::vector<std::complex<double>> JpowerIntegrals::calculate_J_power_integrals(int l, double a, int M, double /*threshold*/) {
    using namespace std::complex_literals;

    std::complex<double> J_0 = 0.25 * std::exp(-0.25i * M_PI) / std::sqrt(M_PI * a) * std::exp(0.25i * static_cast<double>(l * l) / a);

    std::complex<double> beta(0, 0.5 / a);
    auto alpha = static_cast<double>(l) * beta;

    std::vector<std::complex<double>> J = {0.0, J_0};

    for (int m = 0; m < M; m++) {
        std::complex<double> term1 = J[J.size() - 1] * alpha;
        std::complex<double> term2 = J[J.size() - 2] * beta * static_cast<double>(m) / static_cast<double>(m + 2);
        std::complex<double> last  = (term1 + term2) / static_cast<double>(m + 3);
        J.push_back(last);
    }

    J.erase(J.begin());
    return J;
}

void JpowerIntegrals::crop(std::vector<std::complex<double>> &J, double threshold) {
    auto isNegligible = [threshold](const std::complex<double> &c) {
        return std::abs(c.real()) < threshold && std::abs(c.imag()) < threshold;
    };
    J.erase(std::find_if_not(J.rbegin(), J.rend(), isNegligible).base(), J.end());
}

} // namespace mrcpp