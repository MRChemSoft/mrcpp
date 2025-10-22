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

#include "JpowerIntegrals.h"
#include <algorithm> // std::find_if_not
#include <fftw3.h>

namespace mrcpp {

JpowerIntegrals::JpowerIntegrals(double a, int scaling, int M, double threshold) {
    this->scaling = scaling;
    int N = 1 << scaling;
    for (int l = 0; l < N; l++) integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
    for (int l = 1 - N; l < 0; l++) integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
}

/// @brief in progress
/// @param index - interger lying in the interval \f$ [ -2^n + 1, \ldots, 2^n - 1 ] \f$.
/// @return in progress
std::vector<std::complex<double>> &JpowerIntegrals::operator[](int index) {
    if (index < 0) index += integrals.size();
    return integrals[index];
}

std::vector<std::complex<double>> JpowerIntegrals::calculate_J_power_integrals(int l, double a, int M, double threshold) {
    using namespace std::complex_literals;

    std::complex<double> J_0 = 0.25 * std::exp(-0.25i * M_PI) / std::sqrt(M_PI * a) * std::exp(0.25i * static_cast<double>(l * l) / a);
    std::complex<double> beta(0, 0.5 / a);
    auto alpha = static_cast<double>(l) * beta;

    std::vector<std::complex<double>> J = {0.0, J_0};

    for (int m = 0; m < M; m++) {
        std::complex<double> term1 = J[J.size() - 1] * alpha;
        std::complex<double> term2 = J[J.size() - 2] * beta * static_cast<double>(m) / static_cast<double>(m + 2);
        std::complex<double> last = (term1 + term2) / static_cast<double>(m + 3);
        J.push_back(last);
    }

    J.erase(J.begin());
    return J;
}

/// @details Removes negligible elements in \b J until it reaches a considerable value.
void JpowerIntegrals::crop(std::vector<std::complex<double>> &J, double threshold) {
    // Lambda function to check if an element is negligible
    auto isNegligible = [threshold](const std::complex<double> &c) { return std::abs(c.real()) < threshold && std::abs(c.imag()) < threshold; };
    // Remove negligible elements from the end of the vector
    J.erase(std::find_if_not(J.rbegin(), J.rend(), isNegligible).base(), J.end());
}


std::vector<std::complex<double>> DerivativePowerIntegrals::calculate_J_power_integrals(int l, double cut_off, int M, double threshold) {
    using namespace std::complex_literals;
    
    const int N = 1 << this->scaling;
    const int total_N = 2 * N + 1;
    const double two_pi = 2.0 * M_PI;

     // --- 2. Frequency grid (same as np.fft.fftfreq) ---
    std::vector<double> xi_freq(total_N);
    for (int k = 0; k < total_N; ++k) {
        double freq = (k <= total_N / 2) ? k : k - total_N;
        xi_freq[k] = two_pi * freq / total_N;
    }

    // --- 3. Smooth cutoff function chi_cut ---
    std::vector<double> chi(total_N);
    for (int k = 0; k < total_N; ++k) {
        double abs_xi = std::abs(N * xi_freq[k]);
        chi[k] = (abs_xi > cut_off)
            ? std::exp(-std::pow(abs_xi - cut_off, 2))
            : 1.0;
    }

    // --- 4. Initialize FFTW arrays ---
    std::vector<std::complex<double>> f_values(total_N);
    for (int k = 0; k < total_N; ++k)
        f_values[k] = chi[k];  // start with chi_cut(ξ)

    std::vector<std::complex<double>> ifft_output(total_N);
    fftw_plan plan = fftw_plan_dft_1d(
        total_N,
        reinterpret_cast<fftw_complex*>(f_values.data()),
        reinterpret_cast<fftw_complex*>(ifft_output.data()),
        FFTW_BACKWARD,  // inverse FFT
        FFTW_ESTIMATE
    );

    // --- 5. Compute power integrals ---
    std::vector<std::vector<double>> power_integrals;
    power_integrals.emplace_back(total_N, 0.0); // dummy for index 0

    for (int m = 1; m < M; ++m) {
        for (int k = 0; k < total_N; ++k)
            f_values[k] *= std::complex<double>(0.0, xi_freq[k]) / double(m + 1);

        fftw_execute(plan);

        // normalize and store real part
        std::vector<double> integral(total_N);
        for (int k = 0; k < total_N; ++k)
            integral[k] = ifft_output[k].real() / total_N;

        power_integrals.push_back(std::move(integral));
    }

    fftw_destroy_plan(plan);



    //std::complex<double> J_0 = 0.25 * std::exp(-0.25i * M_PI) / std::sqrt(M_PI * a) * std::exp(0.25i * static_cast<double>(l * l) / a);
    //std::complex<double> beta(0, 0.5 / a);
    //auto alpha = static_cast<double>(l) * beta;

    std::vector<std::complex<double>> J = {0.0};
/*
    for (int m = 0; m < M; m++) {
        std::complex<double> term1 = J[J.size() - 1] * alpha;
        std::complex<double> term2 = J[J.size() - 2] * beta * static_cast<double>(m) / static_cast<double>(m + 2);
        std::complex<double> last = (term1 + term2) / static_cast<double>(m + 3);
        J.push_back(last);
    }

    J.erase(J.begin());
*/
    return J;
}

} // namespace mrcpp
