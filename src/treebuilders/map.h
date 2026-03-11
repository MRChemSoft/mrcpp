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

/**
 * @file
 * @brief Nonlinear mapping utilities for multiresolution trees.
 *
 * @details
 * Declares an adaptive routine that applies a user-supplied scalar mapping
 * to a multiresolution function represented by a @ref mrcpp::FunctionTree.
 * The routine produces an output tree whose grid is refined as needed to meet
 * a requested precision.
 *
 * ### What “map” does
 * Given an input scalar field \( f(\mathbf{r}) \) encoded by `inp`, and a
 * scalar-to-scalar function `fmap : ℝ → ℝ`, this routine builds (or refines)
 * the topology of `out` and computes coefficients so that
 * \f[
 *    g(\mathbf{r}) = \mathrm{fmap}\big(f(\mathbf{r})\big)
 * \f]
 * is represented to within the requested tolerance.
 *
 * The mapping is *pointwise* in value space (nonlinear allowed) and the grid
 * refinement is *adaptive*: nodes are split where approximation error indicates
 * additional resolution is required.
 *
 * ### Typical uses
 * - Envelope shaping (e.g., clamp, softplus, \f$x^p\f$).
 * - Nonlinearities inside iterative solvers.
 * - Post-processing fields (e.g., magnitude, thresholding).
 *
 * @note Only the **real** scalar case (`double` coefficients) is declared here.
 * Complex-valued mappings typically require splitting real/imag components
 * explicitly or using dedicated complex routines elsewhere in the library.
 */

#include "trees/FunctionTreeVector.h"

namespace mrcpp {

template <int D, typename T> class FunctionTree;

/**
 * @brief Apply a scalar mapping to a function tree with adaptive refinement.
 *
 * @tparam D Spatial dimension (1–3 typical).
 *
 * @param[in]  prec
 *   Target precision threshold used to control adaptive refinement.
 *   See @p absPrec for interpretation.
 * @param[out] out
 *   Destination tree. Its topology will be enlarged/refined as needed and its
 *   coefficients overwritten with the mapped result.
 *   The tree must be associated with a valid MRA compatible with @p inp.
 * @param[in]  inp
 *   Source tree that encodes the input function \( f(\mathbf{r}) \).
 *   Logically read-only (will not be modified by a correct implementation).
 * @param[in]  fmap
 *   Scalar mapping functor of type `FMap<double,double>` (typically equivalent
 *   to `std::function<double(double)>` or any callable with signature
 *   `double(double)`). It is applied pointwise to sample values of `inp`.
 * @param[in]  maxIter
 *   Maximum number of refinement passes. A negative value (default) requests
 *   unbounded passes until the internal convergence criterion is satisfied
 *   (e.g., no new nodes created or estimated error below @p prec everywhere).
 * @param[in]  absPrec
 *   If `true`, interpret @p prec as an **absolute** tolerance on the local
 *   error indicator. If `false` (default), use a **relative** tolerance,
 *   typically scaled by an estimate of \f$\|f\|\f$ (implementation-defined).
 *
 * @pre
 * - `out` and `inp` share compatible MRAs (same domain, basis order, etc.).
 * - `fmap` must be pure (side-effect free) and thread-safe.
 *
 * @post
 * - `out`’s topology and coefficients represent
 *   \( g(\mathbf{r}) = \mathrm{fmap}(f(\mathbf{r})) \) to within the requested
 *   tolerance, subject to the library’s split criterion.
 *
 * @par Precision semantics
 * - *Absolute mode* (`absPrec=true`): the error indicator is compared directly
 *   to @p prec.
 * - *Relative mode* (`absPrec=false`): the indicator is scaled by a norm of
 *   the input (e.g., tree square-norm), so @p prec represents a relative
 *   threshold.
 *
 * @par Parallelization
 * The routine may exploit OpenMP/MPI internally. Supplying a thread-safe
 * `fmap` is required.
 *
 * @par Exception safety
 * Strong guarantee for `inp`. `out` is modified during execution; in case of
 * failure it may be left partially updated.
 *
 * @par Example
 * @code
 * using Tree = mrcpp::FunctionTree<3,double>;
 * Tree fout(mra), fin(mra);
 * // ... build/project fin ...
 *
 * auto square = [](double x){ return x*x; };
 * mrcpp::map<3>(1e-6, fout, fin, square); // fout ≈ (fin)^2
 * @endcode
 */
template <int D>
void map(double prec,
         FunctionTree<D, double> &out,
         FunctionTree<D, double> &inp,
         FMap<double, double> fmap,
         int maxIter = -1,
         bool absPrec = false);

} // namespace mrcpp
