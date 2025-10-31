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

#include "MRCPP/mrcpp_declarations.h"
#include "utils/math_utils.h"

namespace mrcpp {
/**
 * @file
 * @brief Utilities for inspecting and transforming Multiwavelet (MW) trees.
 *
 * @details
 * This header declares helper routines that operate on MRCPP tree structures:
 * - adaptive refinement decisions based on wavelet norms,
 * - creation of per-scale or flat node tables (Hilbert-ordered),
 * - forward and backward multiwavelet transforms between parent/children
 *   scaling coefficients.
 *
 * Unless otherwise stated, functions are **not** thread-safe; synchronize at
 * a higher level if multiple threads may act on the same tree or buffers.
 */
namespace tree_utils {

/**
 * @brief Decide whether a node should be split (refined) based on its wavelet norm.
 *
 * @tparam D Spatial dimension of the MW tree.
 * @tparam T Coefficient type (`double` or `ComplexDouble`).
 * @param node       Node to be tested.
 * @param prec       Target accuracy (relative by default). Non-positive disables splitting.
 * @param split_fac  Scale-dependent factor. If `> MachineZero`, the threshold is
 *                   scaled by \f$2^{-0.5 \cdot \text{split\_fac} \cdot (s+1)}\f$
 *                   where `s` is the node scale; this makes refinement stricter
 *                   at finer scales.
 * @param abs_prec   When `true`, interpret `prec` as an **absolute** tolerance.
 *                   When `false`, use a **relative** tolerance multiplied by
 *                   \f$\|f\|\f$ (square-norm taken from the owning tree).
 *
 * @return `true` if the node’s wavelet norm exceeds the computed threshold and
 *         the node should be refined; `false` otherwise.
 *
 * @details
 * The decision compares \f$\|\mathbf{w}\|\f$ (node wavelet norm) to a threshold:
 * \f[
 *   \tau = \max(2\,\text{MachinePrec},\;
 *               \text{prec} \times (\text{abs\_prec} ? 1 : \|f\|) \times \text{scale\_fac})
 * \f]
 * where \f$\text{scale\_fac}\f$ is determined by `split_fac` as described above.
 * If the owning tree’s square norm is zero and `abs_prec == false`, a fallback
 * of \f$\|f\|=1\f$ is used.
 */
template <int D, typename T>
bool split_check(const MWNode<D, T> &node, double prec, double split_fac, bool abs_prec);

/**
 * @brief Build a flat, Hilbert-ordered table of all non-root nodes in a tree.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param tree   Input MW tree.
 * @param table  Output vector receiving pointers to all internal and leaf nodes
 *               (root depth 0 is skipped). Nodes are traversed in a Hilbert
 *               space-filling curve order; generator nodes are excluded.
 *
 * @details
 * Useful for linear passes (e.g., I/O, diagnostics, custom sweeps) where a
 * contiguous list of nodes is required.
 */
template <int D, typename T>
void make_node_table(MWTree<D, T> &tree, MWNodeVector<D, T> &table);

/**
 * @brief Build per-scale Hilbert-ordered node tables.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param tree    Input MW tree.
 * @param table   Output vector of vectors. Index `d` stores node pointers whose
 *                depth corresponds to `d - tree.getNNegScales()`. Each inner
 *                vector is Hilbert-ordered; generator nodes are excluded.
 *
 * @details
 * This form is convenient for level-wise processing such as multigrid cycles,
 * visualization, or per-scale statistics.
 */
template <int D, typename T>
void make_node_table(MWTree<D, T> &tree, std::vector<MWNodeVector<D, T>> &table);

/**
 * @brief Forward MW transform: build children scaling coefficients from a parent block.
 *
 * @tparam D Spatial dimension (implemented for 1, 2, 3).
 * @tparam T Coefficient type (`double` or `ComplexDouble`).
 * @param tree            Tree providing filter and arity/meta information.
 * @param coeff_in        Pointer to the parent block (size = `kp1^D` entries),
 *                        laid out in standard MRCPP order.
 * @param coeff_out       Pointer to the destination buffer for **children** blocks.
 *                        This routine writes (or accumulates) into `2^D` child
 *                        blocks separated by `stride` elements each.
 * @param readOnlyScaling If `true`, operate as if only scaling components are
 *                        present (skips mixing with wavelets internally).
 * @param stride          Stride, in elements, between consecutive child blocks
 *                        inside `coeff_out`. Must be at least `kp1^D`.
 * @param overwrite       When `true` (default), assign into `coeff_out`.
 *                        When `false`, accumulate (add) into existing values.
 *
 * @pre
 * - `coeff_out` points to sufficient writable storage:
 *   at least `2^D * stride` elements of type `T`.
 * - `coeff_in` points to at least `kp1^D` elements.
 *
 * @post
 * - The `2^D` children scaling blocks are produced in-place in `coeff_out`.
 *
 * @note
 * Complexity is \f$O(2^D \cdot k^{D+1})\f$ for polynomial order `k` (where `kp1 = k+1`).
 * For `D > 3` the routine is not implemented.
 */
template <int D, typename T>
void mw_transform(const MWTree<D, T> &tree,
                  T *coeff_in,
                  T *coeff_out,
                  bool readOnlyScaling,
                  int stride,
                  bool overwrite = true);

// template <int D, typename T> void mw_transform_back(MWTree<D, T> &tree, T *coeff_in, T *coeff_out, int stride);

/**
 * @brief Backward MW transform (3D specialization): build the parent block from children.
 *
 * @tparam T Coefficient type (`double` or `ComplexDouble`).
 * @param tree      Tree providing filter and arity/meta information.
 * @param coeff_in  Pointer to the concatenated **children** blocks (8 blocks in 3D),
 *                  each of size `kp1^3`, separated by `stride` elements.
 * @param coeff_out Pointer to the **parent** block storage (size `kp1^3`).
 * @param stride    Stride, in elements, between consecutive children blocks.
 *
 * @pre
 * - `coeff_in` provides at least `8 * stride` elements.
 * - `coeff_out` provides at least `kp1^3` writable elements.
 *
 * @post
 * - The parent scaling block is reconstructed into `coeff_out`.
 *
 * @note
 * Only the \f$D=3\f$ variant is provided. Use @ref mw_transform for the forward direction.
 */
template <typename T>
void mw_transform_back(MWTree<3, T> &tree, T *coeff_in, T *coeff_out, int stride);

} // namespace tree_utils
} // namespace mrcpp
