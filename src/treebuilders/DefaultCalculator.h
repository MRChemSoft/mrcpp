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
 * @brief Trivial calculator that clears coefficients and norms on each node.
 *
 * @details
 * Declares @ref mrcpp::DefaultCalculator, a minimal
 * @ref TreeCalculator implementation that:
 *  - iterates over a vector of nodes **sequentially** (no OpenMP),
 *  - for each node, clears coefficient/flag state via
 *    @ref MWNode::clearHasCoefs and resets stored norms via
 *    @ref MWNode::clearNorms.
 *
 * This is useful as a baseline or as a final cleanup pass when no numerical
 * operator needs to be applied but tree state must be normalized.
 */

#include "TreeCalculator.h"

namespace mrcpp {

/**
 * @class DefaultCalculator
 * @brief Minimal calculator that performs per-node cleanup.
 *
 * @tparam D Spatial dimension of the tree.
 * @tparam T Scalar coefficient type (`double`, `std::complex<double>`, …).
 *
 * @details
 * The calculator eschews OpenMP parallelism for its node-vector traversal
 * because the work per node is trivial and sequential iteration is typically
 * faster (lower overhead). If parallel traversal is desired, use or derive
 * from an alternative calculator that enables OpenMP in
 * `calcNodeVector`.
 */
template <int D, typename T>
class DefaultCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Process a vector of nodes sequentially.
     *
     * @param nodeVec Container of node pointers to process.
     *
     * @details
     * Calls @ref calcNode on each entry in order. The method deliberately
     * avoids OpenMP to minimize overhead for very small, constant-time work.
     *
     * @complexity Linear in `nodeVec.size()`.
     */
    void calcNodeVector(MWNodeVector<D, T> &nodeVec) override {
        int nNodes = nodeVec.size();
        for (int n = 0; n < nNodes; n++) { calcNode(*nodeVec[n]); }
    }

private:
    /**
     * @brief Clear coefficient presence flags and stored norms for a node.
     *
     * @param node The node whose local state will be reset.
     *
     * @details
     * - @ref MWNode::clearHasCoefs marks that the node no longer has valid
     *   coefficients.
     * - @ref MWNode::clearNorms removes any cached norm values.
     *
     * This does **not** modify topology (no splitting/merging) and does not
     * change coefficient arrays beyond clearing the "has coefs" state.
     */
    void calcNode(MWNode<D, T> &node) override {
        node.clearHasCoefs();
        node.clearNorms();
    }
};

} // namespace mrcpp