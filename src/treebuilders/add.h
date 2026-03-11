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
 * @file add.h
 * @brief Adaptive linear combination of multiresolution (MW) function trees.
 *
 * @details
 * These routines build an output MW function as a weighted sum of one or more
 * input trees on an adaptively refined grid. The refinement loop is driven by a
 * precision target:
 *  - **Relative precision** (default): refine while local wavelet details are
 *    not small compared to the local function norm.
 *  - **Absolute precision** (`absPrec = true`): refine until local details
 *    fall below a fixed absolute threshold.
 *
 * Unless noted otherwise, all input trees and the output must share the same
 * `MultiResolutionAnalysis` (domain, basis, scales). The output grid is
 * extended as needed; it is not cleared automatically.
 */

namespace mrcpp {

/**
 * @brief Adaptive sum of two MW functions with scalar weights.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (e.g., `double`, `ComplexDouble`).
 *
 * @param[in]  prec      Target build precision (relative by default; see @p absPrec).
 * @param[out] out       Output tree to construct (its grid is extended as needed).
 * @param[in]  a         Scalar weight multiplying @p tree_a.
 * @param[in]  tree_a    First input function tree.
 * @param[in]  b         Scalar weight multiplying @p tree_b.
 * @param[in]  tree_b    Second input function tree.
 * @param[in]  maxIter   Maximum refinement iterations; negative means unbounded.
 * @param[in]  absPrec   If true, interpret @p prec as absolute; else relative.
 * @param[in]  conjugate If true and @p T is complex, apply complex conjugation
 *                       to the second operand (@p tree_b) during accumulation.
 *
 * @details
 * Builds
 * \f[
 *   \text{out} \leftarrow a\,\text{tree\_a} \;+\;
 *   b\,(\,\text{conjugate}\ ?\ \overline{\text{tree\_b}}:\text{tree\_b}\,),
 * \f]
 * on a grid refined to meet @p prec under the chosen precision policy.
 */
template <int D, typename T>
void add(double prec,
         FunctionTree<D, T> &out,
         T a,
         FunctionTree<D, T> &tree_a,
         T b,
         FunctionTree<D, T> &tree_b,
         int maxIter = -1,
         bool absPrec = false,
         bool conjugate = false);

/**
 * @brief Adaptive linear combination from a vector of (coefficient, tree) pairs.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type.
 *
 * @param[in]  prec      Target build precision (relative by default; see @p absPrec).
 * @param[out] out       Output tree to construct.
 * @param[in]  inp       Vector of pairs \f$(\alpha_k, f_k)\f$ (type `FunctionTreeVector`).
 * @param[in]  maxIter   Maximum refinement iterations; negative means unbounded.
 * @param[in]  absPrec   If true, interpret @p prec as absolute; else relative.
 * @param[in]  conjugate If true and @p T is complex, apply complex conjugation
 *                       to all trees except the first one during accumulation.
 *
 * @details
 * Builds
 * \f[
 *   \text{out} \leftarrow \sum_k \alpha_k\, g_k,
 * \f]
 * where \f$g_k = \overline{f_k}\f$ if @p conjugate is true (and \f$k>0\f$ in the
 * complex case), otherwise \f$g_k = f_k\f$. The grid is refined adaptively
 * to satisfy @p prec.
 */
template <int D, typename T>
void add(double prec,
         FunctionTree<D, T> &out,
         FunctionTreeVector<D, T> &inp,
         int maxIter = -1,
         bool absPrec = false,
         bool conjugate = false);

/**
 * @brief Convenience overload: adaptive sum of a list of trees with unit weights.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type.
 *
 * @param[in]  prec      Target build precision (relative by default; see @p absPrec).
 * @param[out] out       Output tree to construct.
 * @param[in]  inp       List of tree pointers; each term is taken with weight 1.
 * @param[in]  maxIter   Maximum refinement iterations; negative means unbounded.
 * @param[in]  absPrec   If true, interpret @p prec as absolute; else relative.
 * @param[in]  conjugate If true and @p T is complex, apply complex conjugation
 *                       to all trees except the first during accumulation.
 *
 * @details
 * Equivalent to the `FunctionTreeVector` overload with all coefficients set to 1.
 */
template <int D, typename T>
void add(double prec,
         FunctionTree<D, T> &out,
         std::vector<FunctionTree<D, T> *> &inp,
         int maxIter = -1,
         bool absPrec = false,
         bool conjugate = false);

} // namespace mrcpp
