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
 * @brief Projection helpers to expand analytic/representable functions on
 *        multiresolution bases (function trees).
 *
 * @details
 * These overloads build or refine an output @ref FunctionTree (or a vector of
 * trees) so that the supplied function(s) are represented to within a target
 * precision. The projection is adaptive: nodes are split where the estimated
 * local error exceeds the tolerance, and coefficients are (re)computed only
 * where needed.
 *
 * **Precision semantics**
 * - If @p absPrec is `false` (default), @p prec is interpreted as a
 *   *relative* tolerance with respect to a suitable global/aggregate norm of
 *   the function (typical L²-relative stopping criterion).
 * - If @p absPrec is `true`, @p prec is treated as an *absolute* tolerance
 *   for local/node-wise thresholds.
 *
 * **Iteration control**
 * - @p maxIter limits the number of refinement passes. Use `-1` for the
 *   default behavior (iterate until the tolerance is reached or the internal
 *   refiner deems the grid converged).
 *
 * **Preconditions and side effects**
 * - @p out is modified in-place (grid may be refined/coarsened; coefficients
 *   are (re)computed).
 * - The @p RepresentableFunction or callable provided by the user must be
 *   well-defined on the domain of @p out’s @ref MultiResolutionAnalysis.
 *
 * @note Implementations typically perform, per node:
 *   1) evaluate the input function on the node’s quadrature/stencil,
 *   2) compute scaling/wavelet coefficients,
 *   3) estimate local error and decide on further splitting,
 *   4) stop when global/local criteria satisfy @p prec or @p maxIter is hit.
 */

#include "MRCPP/mrcpp_declarations.h"
#include "trees/FunctionTreeVector.h"
#include <functional>

namespace mrcpp {

/**
 * @brief Project a @ref RepresentableFunction onto an output function tree.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type (e.g., `double`, `ComplexDouble`).
 *
 * @param prec     Target tolerance (relative by default, see @p absPrec).
 * @param out      Destination @ref FunctionTree; refined and filled in-place.
 * @param inp      Analytic / representable function to project.
 * @param maxIter  Maximum refinement passes (`-1` = default/unlimited).
 * @param absPrec  If `true`, interpret @p prec as an absolute tolerance.
 *
 * @details
 * Builds an adaptive multiresolution representation of @p inp in @p out.
 * Existing content of @p out may be reused and further refined.
 */
template <int D, typename T = double>
void project(double prec,
             FunctionTree<D, T> &out,
             RepresentableFunction<D, T> &inp,
             int maxIter = -1,
             bool absPrec = false);

/**
 * @brief Project a user-supplied callable onto an output function tree.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type (e.g., `double`, `ComplexDouble`).
 *
 * @param prec     Target tolerance (relative by default, see @p absPrec).
 * @param out      Destination @ref FunctionTree; refined and filled in-place.
 * @param func     Callable (e.g., lambda) mapping @ref Coord to @p T.
 * @param maxIter  Maximum refinement passes (`-1` = default/unlimited).
 * @param absPrec  If `true`, interpret @p prec as an absolute tolerance.
 *
 * @details
 * Equivalent to the @ref RepresentableFunction overload, but accepts any
 * `std::function<T(const Coord<D>&)>` (or compatible lambda) as the source.
 */
template <int D, typename T = double>
void project(double prec,
             FunctionTree<D, T> &out,
             std::function<T(const Coord<D> &r)> func,
             int maxIter = -1,
             bool absPrec = false);

/**
 * @brief Project multiple callables into a vector of function trees.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type (e.g., `double`, `ComplexDouble`).
 *
 * @param prec     Target tolerance (relative by default, see @p absPrec).
 * @param out      Destination @ref FunctionTreeVector; each entry refined /
 *                 filled in-place. It is expected to have the same length as
 *                 @p func (one tree per callable).
 * @param func     Collection of callables, each mapping @ref Coord to @p T.
 * @param maxIter  Maximum refinement passes (`-1` = default/unlimited).
 * @param absPrec  If `true`, interpret @p prec as an absolute tolerance.
 *
 * @details
 * Applies the single-tree callable projection to each element, pairing
 * `out[i]` with `func[i]`. All trees should share a compatible
 * @ref MultiResolutionAnalysis.
 */
template <int D, typename T = double>
void project(double prec,
             FunctionTreeVector<D, T> &out,
             std::vector<std::function<T(const Coord<D> &r)>> func,
             int maxIter = -1,
             bool absPrec = false);

} // namespace mrcpp