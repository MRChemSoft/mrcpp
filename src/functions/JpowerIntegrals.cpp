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


// ===================== DerivativePowerIntegrals (added) =====================
#include <complex>
#include <vector>
#include <cmath>

#ifdef MRCPP_HAVE_FFTW
  #include <fftw3.h>
#endif

namespace mrcpp {

DerivativePowerIntegrals::DerivativePowerIntegrals(double cut_off, int scaling, int M, double /*threshold*/)
    : scaling(scaling)
{
    // Build the table indexed by ell (shift), then by power m.
    // We compute power-integral arrays via spectral method and transpose them.
    integrals = calculate_J_power_integrals(cut_off, M /*,threshold*/);
}

/*
 * The returned structure MUST be:
 *   integrals[ell_index][m]  with:
 *     - ell_index runs over ℓ = 0..N-1, then ℓ = 1-N..-1   (total 2N-1 entries)
 *     - m runs from 0..M-1 but index 0 is a dummy (unused), i.e. first used is m=1
 */
std::vector<std::vector<double>>
DerivativePowerIntegrals::calculate_J_power_integrals(double cut_off, int M, double /*threshold*/) {
    using cplx = std::complex<double>;

    const int N       = 1 << this->scaling;   // 2^n
    const int total_N = 2 * N + 1;            // indices ℓ = -N..N
    const double two_pi = 2.0 * M_PI;

    // Frequency grid like numpy.fft.fftfreq scaled by 2π/total_N
    std::vector<double> xi(total_N);
    for (int k = 0; k < total_N; ++k) {
        const int freq = (k <= total_N/2) ? k : k - total_N; // 0,1,...,N, -(N-1),..., -1
        xi[k] = two_pi * static_cast<double>(freq) / static_cast<double>(total_N);
    }

    // Smooth cutoff χ(ξ): 1 for |N ξ| ≤ cut_off, decays (Gaussian) outside
    std::vector<double> chi(total_N);
    for (int k = 0; k < total_N; ++k) {
        const double abs_xi = std::abs(N * xi[k]);
        chi[k] = (abs_xi > cut_off) ? std::exp(-std::pow(abs_xi - cut_off, 2.0)) : 1.0;
    }

    // power_integrals[m][k] will store inverse-DFT values for each m
    // We keep index 0 as a dummy vector to preserve the "first used is m=1" convention.
    std::vector<std::vector<double>> power_integrals;
    power_integrals.emplace_back(total_N, 0.0); // dummy m=0

#ifdef MRCPP_HAVE_FFTW
    // ---------- FFTW path ----------
    std::vector<cplx> f_values(total_N);
    std::vector<cplx> ifft_out(total_N);

    // start with f(ξ) = χ(ξ)
    for (int k = 0; k < total_N; ++k) f_values[k] = chi[k];

    fftw_plan plan = fftw_plan_dft_1d(
        total_N,
        reinterpret_cast<fftw_complex*>(f_values.data()),
        reinterpret_cast<fftw_complex*>(ifft_out.data()),
        FFTW_BACKWARD, // inverse FFT
        FFTW_ESTIMATE
    );

    for (int m = 1; m < M; ++m) {
        // Multiply by (i ξ)/(m+1) in Fourier side each step (recurrence)
        for (int k = 0; k < total_N; ++k) {
            f_values[k] *= cplx(0.0, xi[k]) / static_cast<double>(m + 1);
        }

        // Inverse FFT (unnormalized) -> divide by total_N to normalize
        fftw_execute(plan);

        std::vector<double> inv(total_N);
        for (int k = 0; k < total_N; ++k) inv[k] = ifft_out[k].real() / static_cast<double>(total_N);
        power_integrals.push_back(std::move(inv));
    }

    fftw_destroy_plan(plan);
#else
    // ---------- Portable O(N^2) inverse DFT fallback ----------
    // Keep identical math/normalization to FFTW branch
    std::vector<cplx> F(total_N);
    for (int k = 0; k < total_N; ++k) F[k] = chi[k];

    for (int m = 1; m < M; ++m) {
        // Multiply by (i ξ)/(m+1)
        for (int k = 0; k < total_N; ++k) F[k] *= cplx(0.0, xi[k]) / static_cast<double>(m + 1);

        // IDFT: f[x] = (1/N) * Σ_k F[k] * exp(+i 2π k x / N); our xi already carries 2π/total_N,
        // so exp(+i xi[k] * x). Here x corresponds to the "shift" index ℓ in [-N..N].
        std::vector<double> inv(total_N, 0.0);
        for (int x = 0; x < total_N; ++x) {
            cplx sum = 0.0;
            for (int k = 0; k < total_N; ++k) {
                sum += F[k] * std::exp(cplx(0.0, xi[k] * static_cast<double>(x)));
            }
            inv[x] = (sum.real() / static_cast<double>(total_N));
        }
        power_integrals.push_back(std::move(inv));
    }
#endif // MRCPP_HAVE_FFTW

    // --------- Transpose: make integrals indexed by ℓ (shift) then by m ---------
    // Your calculators expect integrals.size() == 2N-1 with ordering:
    //   ℓ = 0,1,...,N-1, then ℓ = 1-N,...,-1
    // Each element is a vector over m = 0..M-1 (index 0 kept as dummy).
    std::vector<std::vector<double>> integrals_by_shift;
    integrals_by_shift.reserve(2 * N - 1);

    auto idx_from_l = [N](int l) { return l + N; }; // map ℓ ∈ [-N..N] to [0..2N]

    // positive and zero shifts: ℓ = 0..N-1
    for (int l = 0; l < N; ++l) {
        std::vector<double> series;
        series.reserve(power_integrals.size());
        for (const auto& inv : power_integrals) {
            series.push_back(inv[idx_from_l(l)]);
        }
        integrals_by_shift.push_back(std::move(series));
    }
    // negative shifts in your historical order: ℓ = 1-N..-1
    for (int l = 1 - N; l < 0; ++l) {
        std::vector<double> series;
        series.reserve(power_integrals.size());
        for (const auto& inv : power_integrals) {
            series.push_back(inv[idx_from_l(l)]);
        }
        integrals_by_shift.push_back(std::move(series));
    }

    return integrals_by_shift;
}

std::vector<double>& DerivativePowerIntegrals::operator[](int index) {
    // Match JpowerIntegrals negative-index behavior
    if (index < 0) index += static_cast<int>(integrals.size());
    return integrals[index];
}

} // namespace mrcpp