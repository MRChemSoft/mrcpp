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

/*
 * MRCPP — JpowerIntegrals + (optional) DerivativePowerIntegrals
 */

#include "JpowerIntegrals.h"
#include <algorithm>   // std::find_if_not
#include <cmath>
#include <complex>
#include <vector>

#ifdef MRCPP_HAVE_FFTW
  #include <fftw3.h>
#endif

namespace mrcpp {

/* -------------------- JpowerIntegrals (original) -------------------- */

JpowerIntegrals::JpowerIntegrals(double a, int scaling, int M, double threshold) {
    this->scaling = scaling;
    int N = 1 << scaling;
    integrals.reserve(2 * N);
    for (int l = 0; l < N; ++l) {
        integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
    }
    for (int l = 1 - N; l < 0; ++l) {
        integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
    }
}

std::vector<std::complex<double>>& JpowerIntegrals::operator[](int index) {
    if (index < 0) index += static_cast<int>(integrals.size());
    return integrals[static_cast<size_t>(index)];
}

std::vector<std::complex<double>>
JpowerIntegrals::calculate_J_power_integrals(int l, double a, int M, double /*threshold*/) {
    using namespace std::complex_literals;

    const std::complex<double> J0 =
        0.25 * std::exp(-0.25i * M_PI) / std::sqrt(M_PI * a)
      * std::exp(0.25i * static_cast<double>(l * l) / a);

    const std::complex<double> beta(0.0, 0.5 / a);
    const std::complex<double> alpha = static_cast<double>(l) * beta;

    std::vector<std::complex<double>> J;
    J.reserve(static_cast<size_t>(M + 2));
    J.push_back(0.0);
    J.push_back(J0);

    for (int m = 0; m < M; ++m) {
        const std::complex<double> term1 = J.back() * alpha;
        const std::complex<double> term2 = J[J.size() - 2] * beta
                                         * static_cast<double>(m) / static_cast<double>(m + 2);
        const std::complex<double> next  = (term1 + term2) / static_cast<double>(m + 3);
        J.push_back(next);
    }

    J.erase(J.begin()); // drop the leading 0
    return J;
}

void JpowerIntegrals::crop(std::vector<std::complex<double>> &J, double threshold) {
    const auto isNegligible = [threshold](const std::complex<double> &c) {
        return std::abs(c.real()) < threshold && std::abs(c.imag()) < threshold;
    };
    J.erase(std::find_if_not(J.rbegin(), J.rend(), isNegligible).base(), J.end());
}

// -----------------------------------------------------------------------------
// DerivativePowerIntegrals implementation
// -----------------------------------------------------------------------------

DerivativePowerIntegrals::DerivativePowerIntegrals(double cut_off, int scaling, int M, double threshold)
    : scaling(scaling)
{
    // match header: calculate_J_power_integrals(double cut_off, int M, double threshold)
    integrals = calculate_J_power_integrals(cut_off, M, threshold);
}

#ifdef MRCPP_HAVE_FFTW
std::vector<std::vector<double>>
DerivativePowerIntegrals::calculate_J_power_integrals(double cut_off, int M, double /*threshold*/)
{
    const int N       = 1 << this->scaling;       // 2^n
    const int total_N = 2 * N + 1;
    const double two_pi = 2.0 * M_PI;

    // Frequency grid (like numpy.fft.fftfreq, scaled by 2π/total_N)
    std::vector<double> xi_freq(total_N);
    for (int kk = 0; kk < total_N; ++kk) {
        const double base = (kk <= total_N / 2) ? kk : kk - total_N;
        xi_freq[kk] = two_pi * base / static_cast<double>(total_N);
    }

    // Smooth cutoff χ(ξ)
    std::vector<double> chi(total_N);
    for (int kk = 0; kk < total_N; ++kk) {
        const double abs_xi = std::abs(N * xi_freq[kk]);
        chi[kk] = (abs_xi > cut_off)
                  ? std::exp(-std::pow(abs_xi - cut_off, 2.0))
                  : 1.0;
    }

    // FFT buffers
    std::vector<std::complex<double>> f_values(total_N);
    for (int kk = 0; kk < total_N; ++kk) f_values[kk] = chi[kk];

    std::vector<std::complex<double>> ifft_out(total_N);
    fftw_plan plan = fftw_plan_dft_1d(
        total_N,
        reinterpret_cast<fftw_complex*>(f_values.data()),
        reinterpret_cast<fftw_complex*>(ifft_out.data()),
        FFTW_BACKWARD,     // inverse FFT
        FFTW_ESTIMATE
    );

    // Accumulate “power integrals”
    std::vector<std::vector<double>> power_integrals;
    power_integrals.emplace_back(total_N, 0.0);   // dummy index 0 (so results align with m = 1..M-1 below)

    for (int m = 1; m < M; ++m) {
        // multiply by (i * ξ)/(m+1) in frequency space
        for (int kk = 0; kk < total_N; ++kk) {
            f_values[kk] *= std::complex<double>(0.0, xi_freq[kk]) / static_cast<double>(m + 1);
        }

        fftw_execute(plan);

        // normalize inverse FFT and store real part
        std::vector<double> inv(total_N);
        for (int kk = 0; kk < total_N; ++kk) {
            inv[kk] = ifft_out[kk].real() / static_cast<double>(total_N);
        }
        power_integrals.push_back(std::move(inv));
    }

    fftw_destroy_plan(plan);
    return power_integrals;
}
#else
// Fallback when FFTW is not available: correct shape, zeros
std::vector<std::vector<double>>
DerivativePowerIntegrals::calculate_J_power_integrals(double /*cut_off*/, int M, double /*threshold*/)
{
    const int N       = 1 << this->scaling;
    const int total_N = 2 * N + 1;

    std::vector<std::vector<double>> power_integrals;
    power_integrals.emplace_back(total_N, 0.0);   // dummy index 0
    for (int m = 1; m < M; ++m)
        power_integrals.emplace_back(total_N, 0.0);

    return power_integrals;
}
#endif

std::vector<double>& DerivativePowerIntegrals::operator[](int index) {
    if (index < 0) index += static_cast<int>(integrals.size());
    return integrals[index];
}

} // namespace mrcpp