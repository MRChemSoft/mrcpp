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
#include <algorithm>
#include <complex>
#include <cmath>
#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

namespace mrcpp {

/* ------------------------------ JpowerIntegrals ------------------------------ */

JpowerIntegrals::JpowerIntegrals(double a, int scaling, int M, double threshold) {
    this->scaling = scaling;
    int N = 1 << scaling;
    for (int l = 0; l < N; ++l)
        integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
    for (int l = 1 - N; l < 0; ++l)
        integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
}

std::vector<std::complex<double>> &JpowerIntegrals::operator[](int index) {
    if (index < 0) index += integrals.size();
    return integrals[index];
}

std::vector<std::complex<double>>
JpowerIntegrals::calculate_J_power_integrals(int l, double a, int M, double /*threshold*/) {
    using namespace std::complex_literals;

    std::complex<double> J0 = 0.25 * std::exp(-0.25i * M_PI) / std::sqrt(M_PI * a)
                              * std::exp(0.25i * static_cast<double>(l * l) / a);
    std::complex<double> beta(0.0, 0.5 / a);
    std::complex<double> alpha = static_cast<double>(l) * beta;

    std::vector<std::complex<double>> J = {0.0, J0};

    for (int m = 0; m < M; ++m) {
        std::complex<double> term1 = J[J.size() - 1] * alpha;
        std::complex<double> term2 = J[J.size() - 2] * beta * static_cast<double>(m) / static_cast<double>(m + 2);
        std::complex<double> next  = (term1 + term2) / static_cast<double>(m + 3);
        J.push_back(next);
    }

    // Drop the leading 0 so that index 0 corresponds to J0
    J.erase(J.begin());
    return J;
}

void JpowerIntegrals::crop(std::vector<std::complex<double>> &J, double threshold) {
    auto isNeg = [threshold](const std::complex<double> &c) {
        return std::abs(c.real()) < threshold && std::abs(c.imag()) < threshold;
    };
    J.erase(std::find_if_not(J.rbegin(), J.rend(), isNeg).base(), J.end());
}

/* --------------------------- DerivativePowerIntegrals --------------------------- */

DerivativePowerIntegrals::DerivativePowerIntegrals(double cut_off, int scaling, int M, double threshold)
    : scaling(scaling) {
    integrals = calculate_J_power_integrals(cut_off, M, threshold);
}

// Returns integrals indexed by spatial shift l ∈ {-(N-1),..., -1, 0, 1, ..., (N-1)}.
// Each entry holds a vector<double> over m = 0..M (m=0 left as 0.0).
std::vector<std::vector<double>>
DerivativePowerIntegrals::calculate_J_power_integrals(double cut_off, int M, double /*threshold*/) {
    const int N        = 1 << this->scaling;
    const int total_N  = 2 * N + 1;                // FFT grid size
    const double two_pi = 2.0 * M_PI;

    // Frequency grid (radians), like numpy.fft.fftfreq mapped to radians
    std::vector<double> xi(total_N);
    for (int k = 0; k < total_N; ++k) {
        double freq = (k <= total_N / 2) ? k : k - total_N;  // 0..N, -(N-1)..-1
        xi[k] = two_pi * freq / total_N;
    }

    // Smooth cutoff χ(ξ)
    std::vector<double> chi(total_N);
    for (int k = 0; k < total_N; ++k) {
        double ax = std::abs(N * xi[k]);
        chi[k] = (ax > cut_off) ? std::exp(-std::pow(ax - cut_off, 2)) : 1.0;
    }

    // FFT buffers
    std::vector<std::complex<double>> f(total_N), out(total_N);
    for (int k = 0; k < total_N; ++k) f[k] = chi[k];

    fftw_plan plan = fftw_plan_dft_1d(
        total_N,
        reinterpret_cast<fftw_complex*>(f.data()),
        reinterpret_cast<fftw_complex*>(out.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    // Storage per shift l: Lspan = 2N - 1  (shifts from -(N-1) to (N-1))
    const int Lspan = 2 * N - 1;
    std::vector<std::vector<double>> by_l(Lspan, std::vector<double>(M + 1, 0.0));

    auto fft_index_of_l = [total_N](int l) -> int {
        return (l >= 0) ? l : (total_N + l);
    };
    auto slot_of_l = [Lspan, N](int l) -> int {
        return (l >= 0) ? l : (Lspan + l);
    };

    // m = 0 stays zero; fill m = 1..M
    for (int m = 1; m <= M; ++m) {
        for (int k = 0; k < total_N; ++k)
            f[k] *= std::complex<double>(0.0, xi[k]) / double(m + 1);

        fftw_execute(plan);

        // scatter to per-l bins
        for (int l = -(N - 1); l <= (N - 1); ++l) {
            int fft_k = fft_index_of_l(l);
            int slot  = slot_of_l(l);
            by_l[slot][m] = out[fft_k].real() / total_N;
        }
    }

    fftw_destroy_plan(plan);
    return by_l;
}

std::vector<double> &DerivativePowerIntegrals::operator[](int index) {
    if (index < 0) index += integrals.size();
    return integrals[index];
}

} // namespace mrcpp