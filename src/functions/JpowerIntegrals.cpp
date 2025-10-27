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
#include <algorithm> // std::find_if_not

namespace mrcpp {

/**
 * # Class purpose
 * Computes, stores, and provides indexed access to the sequence of
 * “power integrals” \( J_m(l) \) for a range of integer shifts \( l \).
 *
 * In this implementation each sequence \(\{J_m(l)\}_{m=0}^{M}\) is produced by
 * a **three–term recurrence** seeded by a closed form for \(J_0(l)\):
 *
 *   - Seed:
 *     \f[
 *       J_0(l)
 *       = \tfrac{1}{4}\,e^{-i\pi/4}\,\frac{1}{\sqrt{\pi a}}\,
 *         \exp\!\Big( \tfrac{i\,l^2}{4a} \Big)
 *     \f]
 *   - Parameters:
 *     \f[
 *       \beta = \tfrac{i}{2a}, \qquad \alpha = l\,\beta
 *     \f]
 *   - Recurrence (implemented below):
 *     \f[
 *       J_{m+2}
 *       = \frac{\alpha\,J_{m+1} + \frac{m}{m+2}\,\beta\,J_{m}}{m+3},
 *       \qquad m=0,1,\dots
 *     \f]
 *
 * The class builds these sequences for all integer \( l \) in the
 * symmetric range \([-(2^n-1), \dots, -1, 0, \dots, 2^n-1]\),
 * where `n = scaling` and `N = 2^n`. Internally, results are stored as
 * `std::vector<std::complex<double>>` (one vector per shift \(l\)).
 *
 * ## Parameters
 * - `a`         : real positive parameter in the Gaussian-like kernel (see seed).
 * - `scaling`   : defines the number of integer shifts as \(N=2^{\text{scaling}}\).
 * - `M`         : the highest power index — sequences contain \(J_0,\dots,J_M\).
 * - `threshold` : magnitude cutoff used by `crop()` to trim negligible tail values.
 *                 (Note: the current constructor does **not** call `crop()`. You
 *                 may call it manually after construction if you want trimming.)
 *
 * ## Indexing
 * The operator `operator[](int index)` accepts the natural range
 * \([-(2^n-1), \dots, 2^n-1]\). Negative indices are transparently
 * mapped to the underlying zero-based container.
 */
JpowerIntegrals::JpowerIntegrals(double a, int scaling, int M, double threshold) {
    this->scaling = scaling;
    int N = 1 << scaling;                 // N = 2^scaling shifts on the positive side (including 0)
    // Store sequences for l = 0,1,...,N-1
    for (int l = 0; l < N; l++) integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
    // And for l = -(N-1),...,-1 (append after the non-negative ones)
    for (int l = 1 - N; l < 0; l++) integrals.push_back(calculate_J_power_integrals(l, a, M, threshold));
}

/**
 * @brief Random–access to the vector of \f$ \{J_m(l)\}_{m=0}^{M} \f$ for a given shift @p index.
 *
 * @param index Integer shift \(l\) in \([-(2^n-1), \dots, 2^n-1]\).
 *              Negative inputs are internally wrapped to the layout used
 *              by the `integrals` storage.
 * @return Reference to the vector `J` containing `[J_0, J_1, ..., J_M]` for that \(l\).
 *
 * @note This is a non-const overload returning a mutable reference; callers
 *       can modify the stored sequence if needed.
 */
std::vector<std::complex<double>> &JpowerIntegrals::operator[](int index) {
    if (index < 0) index += integrals.size(); // wrap negative l to the back half of the container
    return integrals[index];
}

/**
 * @brief Build the sequence \f$ \{J_m(l)\}_{m=0}^{M} \f$ using the closed-form seed and recurrence.
 *
 * @param l         Integer shift parameter.
 * @param a         Positive real parameter from the analytic form.
 * @param M         Highest power index to compute (inclusive).
 * @param threshold Magnitude threshold (currently not used inside this routine).
 * @return Vector of length \f$ M+1 \f$ with entries \f$ [J_0, J_1, \dots, J_M] \f$.
 *
 * Implementation notes:
 * - We store an initial dummy 0 followed by \(J_0\) so that the recurrence
 *   can read the two previous entries uniformly; we erase the dummy before return.
 * - Complex constants:
 *   * `i` is introduced through `std::complex` literals (`std::complex_literals`).
 *   * \f$ \beta = i/(2a) \f$, \f$ \alpha = l \beta \f$.
 * - Numerical behavior:
 *   * The recurrence is simple and linear; for large |m| or extreme `a` you may
 *     see accumulation of round-off; consider `crop()` afterwards if you know
 *     the tail becomes negligible for your use-case.
 */
std::vector<std::complex<double>> JpowerIntegrals::calculate_J_power_integrals(int l, double a, int M, double /*threshold*/) {
    using namespace std::complex_literals;

    // Seed J0(l) = (1/4) e^{-iπ/4} / sqrt(π a) * exp( i l^2 / (4 a) )
    std::complex<double> J_0 = 0.25 * std::exp(-0.25i * M_PI) / std::sqrt(M_PI * a) * std::exp(0.25i * static_cast<double>(l * l) / a);

    // β = i/(2a)  and  α = l β
    std::complex<double> beta(0, 0.5 / a);
    auto alpha = static_cast<double>(l) * beta;

    // Work buffer: prepend a dummy zero so that J.back() = J_m, J[J.size()-2] = J_{m-1}
    // After the loop we drop the dummy, leaving [J_0, J_1, ..., J_M].
    std::vector<std::complex<double>> J = {0.0, J_0};

    // Three-term recurrence:
    // J_{m+2} = (α J_{m+1} + (m/(m+2)) β J_m) / (m+3),  for m = 0..M-1
    for (int m = 0; m < M; m++) {
        std::complex<double> term1 = J[J.size() - 1] * alpha; // α J_{m+1}
        std::complex<double> term2 = J[J.size() - 2] * beta * static_cast<double>(m) / static_cast<double>(m + 2); // (m/(m+2)) β J_m
        std::complex<double> last  = (term1 + term2) / static_cast<double>(m + 3); // divide by (m+3)
        J.push_back(last); // append J_{m+2}
    }

    // Remove the initial dummy zero so the vector starts with J_0
    J.erase(J.begin());
    return J;
}

/**
 * @brief Trim a sequence by removing small-magnitude values from its tail.
 *
 * @param J         The sequence \f$ [J_0, J_1, \dots] \f$ to be cropped in-place.
 * @param threshold Elements with both |real| and |imag| < threshold are considered negligible.
 *
 * Details:
 *  - Traverses from the end until it finds the first element whose real/imag
 *    magnitude is **not** negligible and erases everything past that point.
 *  - Use this to keep only the “significant” prefix \f$ J_0,\dots,J_{m^\*} \f$
 *    if you know the tail rapidly vanishes for your parameters.
 *
 * @warning The constructor does not call this automatically. If you want
 *          trimmed sequences, call `crop(...)` after construction.
 */
void JpowerIntegrals::crop(std::vector<std::complex<double>> &J, double threshold) {
    // Predicate: element is negligible if both real and imaginary parts are below threshold
    auto isNegligible = [threshold](const std::complex<double> &c) {
        return std::abs(c.real()) < threshold && std::abs(c.imag()) < threshold;
    };
    // Erase the trailing run of negligible entries
    J.erase(std::find_if_not(J.rbegin(), J.rend(), isNegligible).base(), J.end());
}

} // namespace mrcpp
