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
 * @brief Adaptor that targets operator-sensitive regions for refinement.
 *
 * @details
 * This adaptor specializes @c WaveletAdaptor in 2D to refine only those nodes
 * that are (i) aligned with a coordinate axis (translation index 0 along @f$x@f$
 * or @f$y@f$) and (ii) exhibit non-zero **wavelet** content. This pattern is
 * useful when applying kernels/operators whose singular support or strongest
 * variation is concentrated along axes (e.g., banded/operator stencils), so
 * we avoid refining benign regions.
 *
 * The refinement trigger is implemented in @ref OperatorAdaptor::splitNode.
 */

#include "WaveletAdaptor.h"

namespace mrcpp {

/**
 * @class OperatorAdaptor
 * @brief Wavelet-driven 2D refinement biased to axis-aligned nodes.
 *
 * @details
 * A node is marked for splitting iff:
 * - **Axis proximity:** its translation index satisfies @c idx[0]==0 or
 *   @c idx[1]==0 (i.e., the node touches the @f$x@f- or @f$y@f-axis at its scale).
 * - **Wavelet energy present:** at least one non-scaling component has a
 *   positive norm. In 2D, component indices are conventionally:
 *   - 0: scaling (S),
 *   - 1..3: wavelet components (W).
 *
 * Combining these filters refines only the regions that are both “near” the
 * axes and relevant to operator action (non-zero wavelet content), keeping the
 * mesh compact elsewhere.
 *
 * @see WaveletAdaptor
 */
class OperatorAdaptor final : public WaveletAdaptor<2> {
public:
    /**
     * @brief Construct an adaptor with precision and depth controls.
     *
     * @param pr Target precision/tolerance forwarded to the base adaptor.
     * @param ms Maximum scale (upper bound on refinement depth) forwarded to the base adaptor.
     * @param ap If @c true, enable the base adaptor's optional parent-aware
     *           behavior (propagation specifics depend on @ref WaveletAdaptor).
     */
    OperatorAdaptor(double pr, int ms, bool ap = false)
            : WaveletAdaptor<2>(pr, ms, ap) {}

protected:
    /**
     * @brief Decide whether a node should be split.
     *
     * @param node The candidate node.
     * @return @c true if the node lies on an axis (either translation index
     *         is zero) **and** has non-zero wavelet component norm; otherwise
     *         @c false.
     *
     * @details
     * - **Axis check:** @c idx[0]==0 || idx[1]==0.
     * - **Wavelet check:** any component @c i in {1,2,3} has
     *   @c node.getComponentNorm(i) > 0.0.
     *
     * Component 0 (scaling) is intentionally ignored in the wavelet check to
     * avoid refining nodes that carry only low-frequency/scaling content.
     */
    bool splitNode(const MWNode<2> &node) const override {
        const auto &idx = node.getNodeIndex();
        bool chkTransl = (idx[0] == 0 || idx[1] == 0);

        bool chkCompNorm = false;
        for (int i = 1; i < 4; i++) {
            if (node.getComponentNorm(i) > 0.0) chkCompNorm = true;
        }

        return chkTransl && chkCompNorm;
    }
};

} // namespace mrcpp