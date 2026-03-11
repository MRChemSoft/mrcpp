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
 * @brief Adaptor that copies data from one or more source trees into a target tree.
 *
 * @details
 * Declares @ref mrcpp::CopyAdaptor, a lightweight @ref TreeAdaptor that
 * drives adaptive traversal/refinement for copy operations.  The adaptor
 * decides whether to split/visit nodes based on a per-dimension band–width
 * window around the regions populated in the source tree(s) and on a
 * max-scale constraint provided at construction.
 */

#include "TreeAdaptor.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {

/**
 * @class CopyAdaptor
 * @brief Adaptor for reproducing (copying) function-tree structure/coefficients.
 *
 * @tparam D Spatial dimensionality (1–3).
 * @tparam T Coefficient scalar type.
 *
 * @details
 * A @ref TreeAdaptor used by generic tree algorithms to:
 *  - restrict refinement to a maximum scale (`ms`), and
 *  - gate splitting to nodes that fall within a per-dimension *band width*
 *    neighborhood around the non-empty region of one or more **source trees**.
 *
 * This enables efficient copying/subsetting of a function tree (or a
 * @ref FunctionTreeVector) into a new tree while avoiding unnecessary
 * refinement outside the area of interest.
 */
template <int D, typename T>
class CopyAdaptor final : public TreeAdaptor<D, T> {
public:
    /**
     * @brief Construct an adaptor using a single source tree.
     *
     * @param t   Source tree to mirror/copy from.
     * @param ms  Maximum scale (depth) allowed for refinement in the target.
     * @param bw  Pointer to an array of length @c D with per-dimension band
     *            half-widths (in node/grid units). Values control how far
     *            from the source support we keep refining; non-positive
     *            entries are treated as zero.
     *
     * @note The adaptor stores an internal vector view that contains @p t.
     */
    CopyAdaptor(FunctionTree<D, T> &t, int ms, int *bw);

    /**
     * @brief Construct an adaptor using multiple source trees.
     *
     * @param t   Collection of source trees whose union of supports guides
     *            refinement/visitation.
     * @param ms  Maximum scale (depth) allowed for refinement in the target.
     * @param bw  Pointer to an array of length @c D with per-dimension band
     *            half-widths (in node/grid units). See the single-tree
     *            constructor for interpretation.
     */
    CopyAdaptor(FunctionTreeVector<D, T> &t, int ms, int *bw);

private:
    /**
     * @brief Per-dimension refinement band half-widths.
     *
     * @details
     * For dimension @c d, only nodes whose index lies within
     * @c bandWidth[d] boxes of the populated region of the source will be
     * considered for splitting. A value of zero limits refinement strictly to
     * the currently populated footprint.
     */
    int bandWidth[D];

    /**
     * @brief Source tree collection used to drive the copy operation.
     *
     * @details
     * When constructed from a single tree, this vector contains exactly one
     * entry referencing that tree; otherwise it aliases the user-provided
     * vector. No ownership transfer takes place.
     */
    FunctionTreeVector<D, T> tree_vec;

    /**
     * @brief Initialize the @ref bandWidth array from a user buffer.
     *
     * @param bw Pointer to an array of length @c D; negative values are clamped
     *           to zero.
     */
    void setBandWidth(int *bw);

    /**
     * @brief Decide whether a node should be split during traversal.
     *
     * @param node Node under consideration in the *target* tree.
     * @return `true` if the node lies within the refinement window and the
     *         max-scale policy permits further subdivision; `false` otherwise.
     *
     * @details
     * The decision combines:
     *  - the maximum allowed scale passed at construction, and
     *  - the per-dimension band width around the union support of the source
     *    tree(s) stored in @ref tree_vec.
     */
    bool splitNode(const MWNode<D, T> &node) const override;
};

} // namespace mrcpp