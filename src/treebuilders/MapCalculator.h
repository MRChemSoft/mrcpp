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
 * @brief Node-wise nonlinear mapping calculator for multiresolution trees.
 *
 * @details
 * This header defines @ref mrcpp::MapCalculator, a concrete
 * @ref mrcpp::TreeCalculator that applies a *pointwise* value mapping
 * (nonlinear allowed) to the coefficients of @ref mrcpp::FunctionTree nodes.
 *
 * The calculator:
 * 1. Locates the corresponding input node (creating a *temporary copy* when
 *    needed) via `FunctionTree::getNode(idx)`.
 * 2. Brings that input node to **value space**: inverse multi-wavelet (MW)
 *    reconstruction followed by a forward coefficient transform (CV) to obtain
 *    per-point values.
 * 3. Applies the user-supplied functor `fmap : T → T` elementwise.
 * 4. Transforms the output node back to the compressed representation
 *    (CV backward, MW compression), sets flags, and updates norms.
 *
 * The calculator operates node-by-node and is typically orchestrated by the
 * surrounding @ref TreeCalculator driver (which handles traversal, work queues,
 * and adaptive refinement).
 */

#include "TreeCalculator.h"

namespace mrcpp {

/**
 * @class MapCalculator
 * @brief Node-local nonlinear mapping (value transform) on a function tree.
 *
 * @tparam D Spatial dimension (e.g., 1, 2, 3).
 * @tparam T Coefficient/value scalar type (e.g., `double`, `ComplexDouble`).
 *
 * @details
 * `MapCalculator` evaluates a user-provided mapping functor `fmap` on the node
 * samples corresponding to the input tree @p inp and writes the transformed
 * values into the output tree managed by the base @ref TreeCalculator.
 *
 * ### Transform pipeline per node
 * For a given output node `node_o` at index `idx`:
 * - Acquire an **input** node copy: `MWNode<D,T> node_i = func->getNode(idx)`.
 * - Convert coefficients to point samples:
 *   - `node_i.mwTransform(Reconstruction);`
 *   - `node_i.cvTransform(Forward);`
 * - Apply `fmap` to each sample: `coefs_o[j] = fmap(coefs_i[j]);`
 * - Convert back to the compressed representation:
 *   - `node_o.cvTransform(Backward);`
 *   - `node_o.mwTransform(Compression);`
 * - Finalize bookkeeping: `setHasCoefs()` and `calcNorms()`.
 *
 * @note
 * The calculator assumes the **input** and **output** trees are compatible
 * (same MRA order/domain and node layouts as scheduled by the driver). Any
 * missing input nodes are generated on-the-fly by `getNode(idx)` as a *copy*.
 *
 * @warning
 * The mapping functor @p fm **must** be thread-safe and side-effect free,
 * as nodes can be processed in parallel by the base calculator.
 */
template <int D, typename T>
class MapCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Construct a node-mapping calculator.
     *
     * @param fm   Elementwise mapping functor `T → T` (copied/moved in).
     * @param inp  Input function tree providing source coefficients.
     *
     * @pre @p inp is initialized on a valid MRA compatible with the target
     *      output tree managed by the driver.
     */
    MapCalculator(FMap<T, T> fm, FunctionTree<D, T> &inp)
            : func(&inp)
            , fmap(std::move(fm)) {}

private:
    /// Pointer to the input function tree (non-owning).
    FunctionTree<D, T> *func;

    /// Elementwise mapping functor applied to node samples.
    FMap<T, T> fmap;

    /**
     * @brief Compute mapped coefficients for one output node.
     *
     * @param node_o Output node to fill (topology assumed prepared by driver).
     *
     * @details
     * - Fetch input node at the same index (creates a temporary node copy).
     * - Perform MW reconstruction and CV forward transforms on the input copy.
     * - Apply @ref fmap to each sample and write into @p node_o.
     * - Restore compressed representation of @p node_o (CV backward, MW compress).
     * - Mark coefficients present and update node norms.
     *
     * @complexity
     * \f$O(n_\text{coef})\f$ per node (excluding transform costs), where
     * \f$n_\text{coef}\f$ is the number of local coefficients.
     *
     * @thread_safety
     * Independent across nodes when the driver runs in parallel. The functor
     * @ref fmap must be reentrant.
     */
    void calcNode(MWNode<D, T> &node_o) override {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        const int n_coefs = node_o.getNCoefs();
        T *coefs_o = node_o.getCoefs();

        // Obtain a temporary input node copy at the same index.
        MWNode<D, T> node_i = func->getNode(idx);

        // Bring input node to value space (reconstruction → forward CV).
        node_i.mwTransform(Reconstruction);
        node_i.cvTransform(Forward);

        // Apply the non-linear map pointwise.
        const T *coefs_i = node_i.getCoefs();
        for (int j = 0; j < n_coefs; ++j) {
            coefs_o[j] = fmap(coefs_i[j]);
        }

        // Return to compressed representation and finalize bookkeeping.
        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

} // namespace mrcpp