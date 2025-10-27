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

/**
 * @file function_utils.cpp
 *
 * @brief Overlap integrals for (possibly polynomially weighted) Cartesian
 *        Gaussian primitives using the Obara–Saika 1D recurrence.
 *
 * Overview
 * --------
 * - `calc_overlap<D>(GaussFunc<D>, GaussFunc<D>)` computes the D-dimensional
 *   overlap of two separable Cartesian Gaussians by factoring the integral
 *   into a product of 1D overlaps along each Cartesian axis.
 * - The core 1D overlap
 *     \f[
 *       S_{p_a p_b}(x_a,x_b;c_a,c_b) =
 *       \int_{-\infty}^{+\infty}\!(x-x_a)^{p_a}(x-x_b)^{p_b}
 *       e^{-c_a (x-x_a)^2}e^{-c_b (x-x_b)^2}\,dx
 *     \f]
 *   is evaluated by `ObaraSaika_ab`, a compact implementation of
 *   the Obara–Saika recurrence relations.
 *
 * Notation (1D)
 * -------------
 * - Exponents: \f$c_a, c_b > 0\f$.
 * - Powers (angular momenta per axis): \f$p_a, p_b \in \mathbb{N}_0\f$.
 * - Centers: \f$x_a, x_b \in \mathbb{R}\f$.
 * - Composite quantities:
 *     \f[
 *       p = c_a + c_b,\quad
 *       \mu = \frac{c_a c_b}{p},\quad
 *       X_{AB} = x_a - x_b,\quad
 *       X_P = \frac{c_a x_a + c_b x_b}{p},\quad
 *       X_{PA} = X_P - x_a,\quad X_{PB} = X_P - x_b.
 *     \f]
 * - Spherical–spherical overlap seed:
 *     \f[
 *       S_{00} = \sqrt{\frac{\pi}{p}}\;\exp(-\mu X_{AB}^2).
 *     \f]
 *
 * Recurrence (sketch)
 * -------------------
 * Let \f$S_{ij}\f$ denote the overlap with powers \f$(i,j)\f$.
 * The code constructs the first “row” \f$S_{0j}\f$ for \f$j=0..p_b\f$
 * via the \f$X_{PB}\f$ recursion, then generates entries with \f$i>0\f$
 * using relations involving \f$X_{AB}\f$ and \f$X_{PA}\f$.
 * Entries are packed into a 1D array `s_coeff` using a simple linear map.
 *
 * Limits
 * ------
 * - `s_coeff` has fixed size 64; the code comment suggests support up to
 *   combined angular momenta roughly \f$p_a \le 20, p_b \le 20\f$
 *   (so that `power_b + 2*power_a` stays within the array).
 */

#include "function_utils.h"

namespace mrcpp {

// Forward declaration of the 1D core routine (defined below).
namespace function_utils {
double ObaraSaika_ab(int power_a, int power_b,
                     double pos_a, double pos_b,
                     double expo_a, double expo_b);
} // namespace function_utils

/**
 * @brief D-dimensional overlap of two separable Cartesian Gaussians.
 *
 * The D-dimensional overlap factorizes into a product of 1D overlaps
 * along each coordinate axis. Each 1D factor is computed by the
 * Obara–Saika recurrence (`ObaraSaika_ab`).
 *
 * Mathematically:
 * \f[
 *   \langle \mathbf{a} | \mathbf{b} \rangle
 *   =
 *   c_a c_b \prod_{d=1}^{D}
 *   \int_{-\infty}^{+\infty}
 *     (x_d - A_d)^{p_{a,d}}
 *     (x_d - B_d)^{p_{b,d}}
 *     e^{-\alpha_d (x_d - A_d)^2}
 *     e^{-\beta_d  (x_d - B_d)^2}\,dx_d,
 * \f]
 * where `getPower()[d] = p_{*,d}`, `getPos()[d] = A_d or B_d`,
 * `getExp()[d] = α_d or β_d`, and `getCoef()` multiplies at the end.
 *
 * @tparam D Dimensionality (1,2,3,...).
 * @param a First Gaussian primitive (powers, position, exponents, coefficient).
 * @param b Second Gaussian primitive (powers, position, exponents, coefficient).
 * @return Overlap integral value.
 */
template <int D>
double function_utils::calc_overlap(const GaussFunc<D> &a, const GaussFunc<D> &b) {
    double S = 1.0;

    // Multiply 1D overlaps across all Cartesian axes
    for (int d = 0; d < D; d++) {
        S *= ObaraSaika_ab(
            a.getPower()[d], b.getPower()[d],
            a.getPos()[d],   b.getPos()[d],
            a.getExp()[d],   b.getExp()[d]
        );
    }

    // Global prefactor from the two primitives
    S *= a.getCoef() * b.getCoef();
    return S;
}

/**
 * @brief 1D Obara–Saika recurrence for Cartesian Gaussian overlap.
 *
 * Computes
 * \f[
 *   S_{ij} =
 *   \int_{-\infty}^{+\infty}
 *     (x-x_a)^i (x-x_b)^j
 *     e^{-c_a (x-x_a)^2}
 *     e^{-c_b (x-x_b)^2}\,dx,
 * \f]
 * returning the value for \f$i = \texttt{power\_a}\f$ and
 * \f$j = \texttt{power\_b}\f$.
 *
 * Parameters
 * ----------
 * @param power_a  \f$p_a\f$ (non-negative integer power about center @p pos_a)
 * @param power_b  \f$p_b\f$ (non-negative integer power about center @p pos_b)
 * @param pos_a    \f$x_a\f$ (center of the first Gaussian)
 * @param pos_b    \f$x_b\f$ (center of the second Gaussian)
 * @param expo_a   \f$c_a\f$ (exponent of the first Gaussian)
 * @param expo_b   \f$c_b\f$ (exponent of the second Gaussian)
 *
 * Implementation notes
 * --------------------
 * - Forms the composite exponent \f$p=c_a+c_b\f$ and reduced exponent
 *   \f$\mu = c_a c_b / p\f$.
 * - Computes the “product center” \f$X_P = (c_a x_a + c_b x_b)/p\f$
 *   and shift distances \f$X_{PA}=X_P-x_a\f$, \f$X_{PB}=X_P-x_b\f$.
 * - Seeds the recurrence with the spherical–spherical overlap
 *   \f$S_{00} = \sqrt{\pi/p}\,\exp(-\mu (x_a-x_b)^2)\f$.
 * - Builds the first row \f$S_{0j}\f$ for \f$j=0..p_b\f$ using the
 *   forward recurrence in \f$j\f$ (involving \f$X_{PB}\f$ and \f$p\f$).
 * - Extends to \f$i>0\f$ by recurrences that couple \f$S_{i0}\f$, \f$S_{i1}\f$
 *   to previously computed entries and the shifts \f$X_{AB}=x_a-x_b\f$,
 *   \f$X_{PA}\f$.
 *
 * Storage
 * -------
 * - Coefficients are stored in a flat array `s_coeff` with a simple linear
 *   indexing that appends new entries as they are generated:
 *     - indices 0..power_b    : `S_{0,0}, S_{0,1}, ..., S_{0,power_b}`
 *     - then pairs `(S_{1,0}, S_{1,1})`, `(S_{2,0}, S_{2,1})`, ...
 * - The last needed value is at index `power_b + 2*power_a`.
 *
 * @return The requested overlap value `S_{power_a, power_b}`.
 */
double function_utils::ObaraSaika_ab(int power_a, int power_b,
                                     double pos_a, double pos_b,
                                     double expo_a, double expo_b) {
    int i, j;
    double expo_p, mu, pos_p, x_ab, x_pa, x_pb, s_00;

    // Maximum size comment from original author:
    // "The highest angular momentum combination is l=20 for a and b simultaneously"
    // With a flat buffer of length 64, the required index (= power_b + 2*power_a)
    // must be < 64. Keep powers within this bound.
    double s_coeff[64];

    // ---- Composite quantities and seed S_00 ----
    expo_p = expo_a + expo_b;                           // p = c_a + c_b
    mu     = expo_a * expo_b / (expo_a + expo_b);       // μ = c_a c_b / p
    pos_p  = (expo_a * pos_a + expo_b * pos_b) / expo_p;// X_P
    x_ab   = pos_a - pos_b;                             // X_AB
    x_pa   = pos_p - pos_a;                             // X_PA
    x_pb   = pos_p - pos_b;                             // X_PB

    s_00 = pi / expo_p;
    s_00 = std::sqrt(s_00) * std::exp(-mu * x_ab * x_ab); // S_{00}

    // ---- First row: S_{0,j} for j=0..power_b ----
    s_coeff[0] = s_00;               // S_{0,0}
    s_coeff[1] = x_pb * s_00;        // S_{0,1}

    j = 1;
    // Recurrence in j:
    //   S_{0,j+1} = X_PB * S_{0,j} + (j / (2p)) * S_{0,j-1}
    while (j < power_b) {
        s_coeff[j + 1] = x_pb * s_coeff[j] + j * s_coeff[j - 1] / (2.0 * expo_p);
        j++;
    }

    // ---- Bootstrap first two entries with i > 0: S_{1,0}, S_{1,1} ----
    // Relations:
    //   S_{1,0} = S_{0,1} - X_AB * S_{0,0}
    //   S_{1,1} = X_PA * S_{1,0} + (j/(2p)) * S_{0,j}  with j = power_b
    s_coeff[j + 1] = s_coeff[j] - x_ab * s_coeff[j - 1];                         // S_{1,0}
    s_coeff[j + 2] = x_pa * s_coeff[j] + j * s_coeff[j - 1] / (2.0 * expo_p);    // S_{1,1}

    i = 1;
    // ---- General i>0 step: append (S_{i+1,0}, S_{i+1,1}) for i=1..power_a-1 ----
    while (i < power_a) {
        int i_l = j + 2 * i + 1; // index for S_{i+1,0}
        int i_r = j + 2 * i + 2; // index for S_{i+1,1}

        // S_{i+1,0} = S_{i,1} - X_AB * S_{i,0}
        s_coeff[i_l] = s_coeff[i_l - 1] - x_ab * s_coeff[i_l - 2];

        // S_{i+1,1} = X_PA * S_{i,1} + (j * S_{i,0} + i * S_{i-1,0}) / (2p)
        // (the packed indexing below matches these dependencies)
        s_coeff[i_r] = x_pa * s_coeff[i_r - 2] + (j * s_coeff[i_r - 3] + i * s_coeff[i_r - 4]) / (2.0 * expo_p);

        i++;
    }

    // The requested entry is S_{power_a, power_b} at index power_b + 2*power_a.
    return s_coeff[power_b + 2 * power_a];
}

// ---- Explicit template instantiations for common dimensions ----
template double function_utils::calc_overlap<1>(const GaussFunc<1> &a, const GaussFunc<1> &b);
template double function_utils::calc_overlap<2>(const GaussFunc<2> &a, const GaussFunc<2> &b);
template double function_utils::calc_overlap<3>(const GaussFunc<3> &a, const GaussFunc<3> &b);

} // namespace mrcpp