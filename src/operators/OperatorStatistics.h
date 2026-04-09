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

/**
 * @file OperatorStatistics.h
 * @brief Thread-aware counters and summaries for multiwavelet operator application.
 *
 * @details
 * This helper aggregates lightweight statistics collected while applying
 * operators to multiwavelet nodes. For performance and thread-safety, counts
 * are first accumulated in per-thread storage and later merged into global
 * totals using @ref flushNodeCounters().
 *
 * Tracked quantities:
 *  - Total number of destination (*f*) nodes where an operator was applied.
 *  - Total number of source (*g*) nodes evaluated.
 *  - Total number of destination nodes marked as “generalized”.
 *  - An 8×8 histogram of component-pair usages (indexed by `(ft, gt)`).
 *
 * The class intentionally avoids synchronization primitives inside hot loops;
 * callers should invoke @ref flushNodeCounters() at safe points to consolidate
 * results and reset per-thread buffers.
 */

#pragma once

#include <Eigen/Core>
#include <iomanip>

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @class OperatorStatistics
 * @brief Collects and reports counters during operator application.
 *
 * @note
 * - Per-thread counters are sized using @c mrcpp_get_max_threads().
 * - Use the stream operator to print a human-readable summary.
 */
class OperatorStatistics final {
public:
    /// Construct an empty statistics object with per-thread accumulators.
    OperatorStatistics();

    /// Release all dynamically allocated per-thread buffers and histograms.
    ~OperatorStatistics();

    /**
     * @brief Consolidate per-thread counters into global totals and reset locals.
     *
     * @details
     * After calling this, @c totFCount, @c totGCount, @c totGenCount and
     * @c totCompCount reflect all work since the previous flush, and the
     * per-thread buffers are zeroed.
     */
    void flushNodeCounters();

    /**
     * @brief Increment destination (*f*)-node counters for the current thread.
     * @tparam D Spatial dimension of the node.
     * @tparam T Coefficient/value type stored by the node.
     * @param fNode Destination node being updated.
     * @param ft    Destination component bitfield (0–7).
     * @param gt    Source component bitfield (0–7).
     *
     * @details
     * Increments the per-thread f-node count, updates the (ft,gt) entry of the
     * per-thread 8×8 histogram, and increments the generalized-node count if
     * @c fNode.isGenNode() returns true.
     */
    template <int D, typename T>
    void incrementFNodeCounters(const MWNode<D, T> &fNode, int ft, int gt);

    /**
     * @brief Increment source (*g*)-node counters for the current thread.
     * @tparam D Spatial dimension of the node.
     * @tparam T Coefficient/value type stored by the node.
     * @param gNode Source node being processed (unused; for interface symmetry).
     */
    template <int D, typename T>
    void incrementGNodeCounters(const MWNode<D, T> &gNode);

    /// Print a summary of accumulated totals and the component histogram.
    friend std::ostream &operator<<(std::ostream &o, const OperatorStatistics &os) { return os.print(o); }

protected:
    int nThreads;                                ///< Number of worker threads.
    int totFCount;                               ///< Global total of applied *f*-nodes.
    int totGCount;                               ///< Global total of processed *g*-nodes.
    int totGenCount;                             ///< Global total of applied generalized nodes.

    int *fCount;                                 ///< Per-thread *f*-node counters (size = nThreads).
    int *gCount;                                 ///< Per-thread *g*-node counters (size = nThreads).
    int *genCount;                               ///< Per-thread generalized-node counters (size = nThreads).

    Eigen::Matrix<int, 8, 8> *totCompCount;      ///< Global (ft,gt) 8×8 usage histogram.
    Eigen::Matrix<int, 8, 8> **compCount;        ///< Per-thread 8×8 usage histograms.

    /// Internal pretty-printer used by the stream operator.
    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp