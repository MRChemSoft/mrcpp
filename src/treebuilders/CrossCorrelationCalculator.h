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
 * @brief Cross-correlation tree calculator for 2D functions with a 1D kernel.
 *
 * @details
 * Declares @ref mrcpp::CrossCorrelationCalculator, a concrete
 * @ref TreeCalculator that evaluates a (discrete) cross–correlation between a
 * two–dimensional multiresolution function and a one–dimensional kernel,
 * typically along one axis of each 2D node.  The implementation leverages a
 * @ref CrossCorrelationCache to reuse banded operator data and reduce
 * per–node setup costs across the traversal.
 */

#include "TreeCalculator.h"
#include "core/CrossCorrelationCache.h"

namespace mrcpp {

/**
 * @class CrossCorrelationCalculator
 * @brief Applies a cached cross–correlation with a 1D kernel to a 2D tree.
 *
 * @details
 * This calculator specializes @ref TreeCalculator for 2D nodes
 * (`TreeCalculator<2>`).  During traversal, @ref calcNode pulls the relevant
 * coefficient band(s) from the current node, applies a cross–correlation with
 * the supplied 1D @ref kernel, and writes the result back to the destination
 * tree/state managed by the base calculator.
 *
 * ### Design notes
 * - The kernel is provided as a `FunctionTree<1>` and is **not owned** by the
 *   calculator (the caller must guarantee its lifetime).
 * - Internally, @ref applyCcc parametrizes on the scalar coefficient type
 *   (`double`, `std::complex<double>`, …) via the template parameter `T` of
 *   @ref CrossCorrelationCache, enabling reuse for both real and complex trees.
 * - A @ref CrossCorrelationCache is used to memoize structure- and
 *   bandwidth-dependent intermediates (e.g., band shapes, transforms) so that
 *   repeated applications across many nodes are efficient.
 */
class CrossCorrelationCalculator final : public TreeCalculator<2> {
public:
    /**
     * @brief Construct a calculator using a given 1D kernel.
     *
     * @param k Reference to a 1D function tree representing the correlation
     *          kernel. The pointer is stored; the object must outlive the
     *          calculator.
     *
     * @note No ownership is transferred; `k` must remain valid for the entire
     *       calculation.
     */
    CrossCorrelationCalculator(FunctionTree<1> &k)
            : kernel(&k) {}

private:
    /**
     * @brief Non-owning pointer to the 1D kernel used in the cross–correlation.
     */
    FunctionTree<1> *kernel;

    /**
     * @brief Compute the cross–correlated output for a single 2D node.
     *
     * @param node The node to process within the current output tree. The node's
     *             scale/index determine the coefficient bands to read/write.
     *
     * @details
     * This override fetches the node's input coefficients (from the source tree
     * configured in the base @ref TreeCalculator), prepares/cache-reuses the
     * operator band via a @ref CrossCorrelationCache, and applies the
     * correlation along the appropriate axis. The resulting coefficients are
     * accumulated into the node's output buffer.
     *
     * @warning The method assumes that the base calculator has already
     *          orchestrated any required refinement and that input/output tree
     *          storage is valid for @p node.
     */
    void calcNode(MWNode<2> &node) override;

    /**
     * @brief Apply the cached cross–correlation to a node with a concrete scalar type.
     *
     * @tparam T Coefficient scalar type of the node (e.g., `double`,
     *           `std::complex<double>`).
     * @param node The target node to which the operator is applied.
     * @param ccc  A cross–correlation cache specialized for @p T that provides
     *             band sizes, temporary buffers, and any preassembled operator
     *             pieces required for efficient application.
     *
     * @details
     * Performs the type-specific core computation:
     *  - obtains the relevant coefficient band(s) from @p node,
     *  - uses @p ccc to assemble or retrieve the required operator slice
     *    derived from @ref kernel,
     *  - applies the correlation (with the proper band width and alignment),
     *  - writes/accumulates the result into the node's output coefficients.
     *
     * The separation into this templated helper allows the public
     * @ref calcNode to dispatch based on the underlying node coefficient type
     * without duplicating logic.
     */
    template <int T>
    void applyCcc(MWNode<2> &node, CrossCorrelationCache<T> &ccc);
};

} // namespace mrcpp