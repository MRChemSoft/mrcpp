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

#include "trees/MWNode.h"

namespace mrcpp {

/**
 * @class TreeCalculator
 * @brief Abstract base for per-node computations on @ref MWTree.
 *
 * @tparam D Spatial dimension of the multiwavelet tree.
 * @tparam T Coefficient value type (e.g. `double`, `ComplexDouble`).
 *
 * @details
 * A `TreeCalculator` defines how to **evaluate/update a single node**
 * (via the pure-virtual @ref calcNode) and provides utilities to apply
 * that logic over a set of nodes, possibly in parallel.
 *
 * Typical usage (in conjunction with @ref TreeBuilder):
 *  - derive a calculator and implement @ref calcNode,
 *  - obtain a worklist with @ref getInitialWorkVector,
 *  - call @ref calcNodeVector to process all nodes,
 *  - optionally override @ref postProcess for statistics/timers.
 *
 * ### Parallelism
 * @ref calcNodeVector uses OpenMP (if enabled) with `schedule(guided)`
 * and a thread count provided by the `mrcpp_get_num_threads()` macro.
 * Implementations of @ref calcNode must be **thread-safe** w.r.t. other
 * nodes in the worklist. Avoid shared mutable state unless properly
 * synchronized.
 */
template <int D, typename T>
class TreeCalculator {
public:
    /// @brief Default constructor.
    TreeCalculator() = default;

    /// @brief Virtual destructor.
    virtual ~TreeCalculator() = default;

    /**
     * @brief Build the initial list of nodes to process.
     *
     * @param[in,out] tree The tree whose nodes should be evaluated.
     * @return Heap-allocated vector of node pointers representing the initial
     *         work set (typically the **current leaf nodes**).
     *
     * @details
     * The default implementation returns a copy of the tree's end-node table
     * (`tree.copyEndNodeTable()`). Callers are responsible for deleting the
     * returned container when done.
     */
    virtual MWNodeVector<D, T>* getInitialWorkVector(MWTree<D, T> &tree) const {
        return tree.copyEndNodeTable();
    }

    /**
     * @brief Evaluate all nodes in @p nodeVec (parallelized when available).
     *
     * @param[in,out] nodeVec Container of node pointers to be processed.
     *
     * @details
     * Invokes @ref calcNode for each entry. Uses OpenMP with guided scheduling
     * and `mrcpp_get_num_threads()` to determine the thread count.
     * After processing all nodes, calls @ref postProcess once.
     *
     * @note The container is treated as read-only regarding its topology;
     *       implementations of @ref calcNode should not insert/remove nodes.
     */
    virtual void calcNodeVector(MWNodeVector<D, T> &nodeVec) {
#pragma omp parallel shared(nodeVec) num_threads(mrcpp_get_num_threads())
        {
            int nNodes = nodeVec.size();
#pragma omp for schedule(guided)
            for (int n = 0; n < nNodes; n++) {
                MWNode<D, T> &node = *nodeVec[n];
                calcNode(node);
            }
        }
        postProcess();
    }

protected:
    /**
     * @brief Perform the calculator's core work on a single node.
     *
     * @param[in,out] node Target node. Implementations typically:
     *  - ensure transforms are in the correct space (MW/CV) as needed,
     *  - compute/update coefficients and derived norms/flags,
     *  - leave the node in a consistent state for subsequent passes.
     *
     * @warning This method is called concurrently on different nodes.
     *          Do not mutate shared global state without synchronization.
     */
    virtual void calcNode(MWNode<D, T> &node) = 0;

    /**
     * @brief Optional hook executed once after @ref calcNodeVector finishes.
     *
     * @details
     * Override to flush accumulators, update statistics, print timers, etc.
     * Default implementation is a no-op.
     */
    virtual void postProcess() {}
};

} // namespace mrcpp