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

#pragma once

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

namespace mrcpp {

/** @class JpowerIntegrals
 *
 * @brief Precompute and store families of power–type integrals \f$ \{\widetilde J_m(l,a)\}_{m\ge 0} \f$
 * for integer shifts \f$ l \f$, used in the Schrödinger time–evolution operator.
 *
 * @details
 * This helper class generates the sequences
 * \f$ \big(\widetilde J_m(l,a)\big)_{m=0}^{M} \f$ for a finite set of integer
 * translations \f$ l \in \{-(2^n-1),\ldots,-1,0,1,\ldots,2^n-1\} \f$, where
 * \f$ n=\texttt{scaling} \f$ and \f$ a>0 \f$ is the time–scaled parameter
 * (typically \f$ a = t\,\mathfrak N^2 = t\,4^{\mathfrak n} \f$).
 *
 * The integrals appear in the expansion of the (matrix–valued) operator
 * \f[
 *   \big[ \sigma_l^{\mathfrak n} \big]_{pj}(a)
 *   =
 *   \sum_{k=0}^{\infty} C_{jp}^{2k}\,
 *   \widetilde J_{\,2k + j + p}(l,a),
 * \f]
 * where the scalar building blocks are
 * \f[
 *   \widetilde J_m(l,a)
 *   =
 *   \frac{e^{i\frac{\pi}{4}(m-1)}}{2\pi\,(m+2)!}
 *   \int_{\mathbb R}
 *     \exp\!\Big(
 *       \rho\,l\,e^{i\pi/4} - a\,\rho^2
 *     \Big)\,
 *     \rho^m\, d\rho .
 * \f]
 *
 * In the code, \f$ \widetilde J_m \f$ are produced by the three–term recurrence
 * (valid for \f$ m=0,1,2,\ldots \f$)
 * \f[
 *   \widetilde J_{m+1}
 *   =
 *   \frac{i}{2a\,(m+3)}\left(
 *     l\,\widetilde J_m + \frac{m}{m+2}\,\widetilde J_{m-1}
 *   \right),
 *   \qquad \widetilde J_{-1}=0,
 * \f]
 * with the closed–form seed
 * \f[
 *   \widetilde J_0(l,a)
 *   =
 *   \frac{e^{-i\pi/4}}{4\sqrt{\pi a}}\,
 *   \exp\!\left(\frac{i\,l^2}{4a}\right).
 * \f]
 *
 * ### Storage layout
 * For convenience of iteration, the container `integrals` is filled
 * in the following order:
 * \f[
 *   l = 0, 1, \ldots, 2^n-1,\; 1-2^n, 2-2^n, \ldots, -2, -1.
 * \f]
 * Each entry is a vector
 * \code
 *   integrals[k] == { J_0(l), J_1(l), ..., J_M(l) }
 * \endcode
 * of complex values for a fixed shift \f$ l \f$.
 *
 * ### Intended use
 * - Construct once for a given \f$ a \f$, \f$ n \f$ and \f$ M \f$.
 * - Access the sequence for a particular shift via `operator[](l)`.
 * - Combine with precomputed correlation coefficients \f$ C_{jp}^{2k} \f$
 *   to assemble \f$ [\sigma_l^{\mathfrak n}]_{pj}(a) \f$.
 *
 * @note The class offers an internal @ref crop routine to trim negligible
 * tail entries of a sequence (based on a magnitude threshold). Whether and
 * when cropping is used is an implementation detail; sequences are always
 * returned in full length \f$ M\!+\!1 \f$ from the constructor path.
 */
class JpowerIntegrals
{
public:
    /// @brief Construct and precompute all \f$ \widetilde J_m(l,a) \f$ for
    ///        \f$ l\in[-(2^n-1),\ldots,2^n-1] \f$ and \f$ m=0,\ldots,M \f$.
    ///
    /// @param a         Time–scaled parameter (typically \f$ a=t\,4^{\mathfrak n} \f$), must be positive.
    /// @param scaling   Level \f$ n \f$ that defines \f$ N=2^n \f$ distinct nonnegative shifts
    ///                  (the negative ones are added symmetrically after them).
    /// @param M         Highest power index in the sequence (inclusive). Each stored vector
    ///                  has length \f$ M+1 \f$ starting at \f$ m=0 \f$.
    /// @param threshold Magnitude cutoff used by the private @ref crop routine to remove
    ///                  negligible tail entries (if cropping is applied internally).
    ///
    /// @details
    /// The internal ordering of the outer container is
    /// \f$ l=0,1,\ldots,2^n-1, 1-2^n,\ldots,-1 \f$. This ordering is mirrored by
    /// the @ref operator[] which will map negative indices to the appropriate
    /// position of the storage.
    JpowerIntegrals(double a, int scaling, int M, double threshold = 1.0e-15);
    //JpowerIntegrals(const JpowerIntegrals& other);

    /// @brief Scaling level \f$ n \f$ (kept for reference; not used directly in lookups).
    int scaling;

    /// @brief Container of sequences \f$ \{\widetilde J_m(l,a)\}_{m=0}^M \f$ for all shifts \f$ l \f$.
    /// Each element is a vector of length \f$ M+1 \f$:
    /// \code
    ///   integrals[idx_for_l] = { J_0(l), J_1(l), ..., J_M(l) }
    /// \endcode
    std::vector<std::vector<std::complex<double>>> integrals;

    /// @brief Mutable access to the precomputed sequence for a given shift \f$ l \f$.
    ///
    /// @param index Integer shift \f$ l \in [-(2^n-1), \ldots, 2^n-1] \f$.
    /// @return Reference to the vector \f$ [J_0(l), \ldots, J_M(l)] \f$.
    ///
    /// @details
    /// Negative indices are transparently remapped to the internal storage order
    /// (see the constructor’s documentation). This allows natural use like `obj[-3]`.
    std::vector<std::complex<double>> & operator[](int index);

private:
    /// @brief Build one full sequence \f$ \{\widetilde J_m(l,a)\}_{m=0}^M \f$ for a fixed shift @p l.
    ///
    /// @param l         Shift index.
    /// @param a         Time–scaled parameter (positive).
    /// @param M         Highest power index.
    /// @param threshold Magnitude cutoff passed to @ref crop (if enabled).
    ///
    /// @return A vector with entries \f$ [J_0(l), J_1(l), \ldots, J_M(l)] \f$.
    ///
    /// @details
    /// The routine uses the closed–form seed \f$ \widetilde J_0(l,a) \f$ and the
    /// recurrence relation to fill the sequence up to \f$ m=M \f$.
    std::vector<std::complex<double>> calculate_J_power_integrals(int l, double a, int M, double threshold);

    /// @brief Remove negligible tail entries from a sequence in place.
    ///
    /// @param J         The sequence to be cropped (modified in place).
    /// @param threshold Entries with both real and imaginary parts below @p threshold
    ///                  in absolute value are considered negligible.
    ///
    /// @details
    /// Cropping can be used to shrink \f$ [J_0,\ldots,J_M] \f$ to
    /// \f$ [J_0,\ldots,J_{m^\*}] \f$ once the tail has decayed under the requested tolerance.
    void crop(std::vector<std::complex<double>> & J, double threshold);
};

} // namespace mrcpp