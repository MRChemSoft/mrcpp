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
 * @brief Adaptive refinement policy for multiplying two multiresolution trees.
 *
 * @details
 * This header defines @ref mrcpp::MultiplicationAdaptor, a concrete
 * @ref mrcpp::TreeAdaptor that decides whether an output node should be split
 * (refined) when computing the **pointwise product** of two
 * @ref mrcpp::FunctionTree objects.
 *
 * The decision is based on inexpensive *a priori* estimates derived from each
 * input node's **maximum scaling** and **maximum wavelet** norms.
 * For corresponding nodes at the same index:
 * - let \f$S_i = \sqrt{\text{max scaling square norm of tree } i}\f$,
 * - let \f$W_i = \sqrt{\text{max wavelet square norm of tree } i}\f$,
 *   for \f$i \in \{0,1\}\f$.
 *
 * The product's wavelet contribution is estimated as
 * \f[
 *   \text{multNorm} \approx W_0 S_1 \;+\; W_1 S_0 \;+\; W_0 W_1 .
 * \f]
 * If `multNorm` exceeds the user threshold `prec`, the adaptor requests a split
 * **unless** both input nodes are already leaves, thereby preventing refinement
 * beyond the resolution of the inputs.
 */

#include "TreeAdaptor.h"
#include "trees/FunctionTreeVector.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @class MultiplicationAdaptor
 * @brief Refinement rule for the product of two function trees.
 *
 * @tparam D Spatial dimension (e.g., 1, 2, 3).
 * @tparam T Coefficient/value scalar type (e.g., `double`, `ComplexDouble`).
 *
 * @details
 * The adaptor is typically used to drive construction of an output grid for
 * \f$f_0 \cdot f_1\f$. At each node index it:
 * 1. Retrieves the corresponding nodes from both input trees.
 * 2. Computes the square-rooted maximum scaling and wavelet norms
 *    \f$(S_i, W_i)\f$.
 * 3. Forms the estimate
 *    \f$\text{multNorm} = W_0 S_1 + W_1 S_0 + W_0 W_1\f$.
 * 4. Requests a split if `multNorm > prec` and at least one input node is not
 *    a leaf. This effectively avoids refining deeper than either input grid,
 *    because when both inputs have zero wavelet contribution at a node,
 *    `multNorm == 0` and no further refinement is triggered.
 *
 * @note The input vector @ref trees must contain **exactly two** trees; a
 * runtime error is emitted otherwise.
 *
 * @warning The adaptor reads norms from the input trees during `splitNode`.
 * The member @ref trees is `mutable` to allow this from a `const` context.
 */
template <int D, typename T>
class MultiplicationAdaptor : public TreeAdaptor<D, T> {
public:
    /**
     * @brief Construct the multiplication refinement rule.
     *
     * @param pr  Refinement threshold \f$\text{prec}\f$ for the multNorm estimate.
     * @param ms  Maximum scale/depth hint passed to @ref TreeAdaptor base.
     * @param t   The pair of input trees as a @ref FunctionTreeVector.
     *
     * @pre `t.size() == 2`
     */
    MultiplicationAdaptor(double pr, int ms, FunctionTreeVector<D, T> &t)
            : TreeAdaptor<D, T>(ms)
            , prec(pr)
            , trees(t) {}

    ~MultiplicationAdaptor() override = default;

protected:
    /// Refinement threshold used against `multNorm`.
    double prec;

    /**
     * @brief Input trees used to estimate the product's local complexity.
     *
     * @details
     * Marked `mutable` so that `splitNode` can retrieve node views from a
     * `const` context without implying logical modification.
     */
    mutable FunctionTreeVector<D, T> trees;

    /**
     * @brief Decide whether an output node should be split.
     *
     * @param node The (prospective) output node whose index determines
     *             which input nodes are inspected.
     * @return `true` if `multNorm > prec` **and** at least one of the input
     *         nodes is not a leaf; otherwise `false`.
     *
     * @details
     * - Retrieves the corresponding nodes from both input trees at the same
     *   @ref NodeIndex.
     * - Computes
     *   \f$S_i=\sqrt{\text{max scaling square norm}},\;
     *     W_i=\sqrt{\text{max wavelet square norm}}\f$.
     * - Forms \f$\text{multNorm}=W_0 S_1 + W_1 S_0 + W_0 W_1\f$.
     * - Triggers refinement when the estimate exceeds @ref prec, except when
     *   both inputs are already leaf nodes at this index.
     *
     * @throws Emits `MSG_ERROR` if `trees.size() != 2`.
     */
    bool splitNode(const MWNode<D, T> &node) const override {
        if (this->trees.size() != 2) MSG_ERROR("Invalid tree vec size: " << this->trees.size());

        auto &pNode0 = get_func(trees, 0).getNode(node.getNodeIndex());
        auto &pNode1 = get_func(trees, 1).getNode(node.getNodeIndex());

        // Square roots convert stored square norms to norms.
        double maxW0 = std::sqrt(pNode0.getMaxWSquareNorm());
        double maxW1 = std::sqrt(pNode1.getMaxWSquareNorm());
        double maxS0 = std::sqrt(pNode0.getMaxSquareNorm());
        double maxS1 = std::sqrt(pNode1.getMaxSquareNorm());

        // Estimated wavelet contribution in the product node.
        double multNorm = maxW0 * maxS1 + maxW1 * maxS0 + maxW0 * maxW1;

        // Never refine beyond both input grids' leaf level.
        if (multNorm > this->prec && !(pNode0.isLeafNode() && pNode1.isLeafNode())) {
            return true;
        } else {
            return false;
        }
    }
};

} // namespace mrcpp