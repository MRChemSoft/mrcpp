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
 * @file OperatorStatistics.cpp
 * @brief Implementation of lightweight counters and summaries used during
 *        multiwavelet operator application.
 *
 * @details
 * This module aggregates per-thread counters while applying operators
 * to multiwavelet nodes. It records:
 *  - Number of *g*-nodes (source nodes) computed.
 *  - Number of *f*-nodes (destination nodes) where an operator was applied.
 *  - Number of *generalized* destination nodes (as reported by MWNode::isGenNode()).
 *  - A small 8×8 histogram of applications by component pair (ft, gt),
 *    where `ft` and `gt` are component bitfields.
 *
 * Thread-local storage is used to avoid contention in hot loops; use
 * flushNodeCounters() to accumulate into totals and reset local counters.
 */

#include "OperatorStatistics.h"
#include "trees/MWNode.h"

using namespace Eigen;

namespace mrcpp {

/**
 * @brief Construct an empty statistics object with per-thread accumulators.
 *
 * @details
 * Allocates:
 *  - @c totCompCount: global 8×8 histogram (zero-initialized).
 *  - Per-thread scalar counters (@c fCount, @c gCount, @c genCount).
 *  - Per-thread 8×8 component histograms (@c compCount[i]).
 *
 * The number of threads is discovered via mrcpp_get_max_threads().
 */
OperatorStatistics::OperatorStatistics()
        : nThreads(mrcpp_get_max_threads())
        , totFCount(0)
        , totGCount(0)
        , totGenCount(0)
        , fCount(nullptr)
        , gCount(nullptr)
        , genCount(nullptr)
        , totCompCount(nullptr)
        , compCount(nullptr) {

    this->totCompCount = new Matrix<int, 8, 8>;
    this->totCompCount->setZero();

    this->fCount = new int[this->nThreads];
    this->gCount = new int[this->nThreads];
    this->genCount = new int[this->nThreads];
    this->compCount = new Matrix<int, 8, 8> *[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->compCount[i] = new Matrix<int, 8, 8>;
        this->compCount[i]->setZero();
        this->fCount[i] = 0;
        this->gCount[i] = 0;
        this->genCount[i] = 0;
    }
}

/**
 * @brief Destroy statistics and free all dynamically allocated arrays.
 */
OperatorStatistics::~OperatorStatistics() {
    for (int i = 0; i < this->nThreads; i++) { delete this->compCount[i]; }
    delete[] this->compCount;
    delete[] this->fCount;
    delete[] this->gCount;
    delete[] this->genCount;
    delete totCompCount;
}

/**
 * @brief Accumulate all per-thread counters into totals and reset locals.
 *
 * @details
 * After this call:
 *  - @c totFCount, @c totGCount, and @c totGenCount are increased by the
 *    sums over all threads.
 *  - @c totCompCount is incremented by each thread-local 8×8 histogram.
 *  - All per-thread counters/histograms are reset to zero.
 */
void OperatorStatistics::flushNodeCounters() {
    for (int i = 0; i < this->nThreads; i++) {
        this->totFCount += this->fCount[i];
        this->totGCount += this->gCount[i];
        this->totGenCount += this->genCount[i];
        *this->totCompCount += *this->compCount[i];
        this->fCount[i] = 0;
        this->gCount[i] = 0;
        this->genCount[i] = 0;
        this->compCount[i]->setZero();
    }
}

/**
 * @brief Increment the *g*-node usage counter for the current thread.
 *
 * @tparam D Spatial dimension of the node.
 * @tparam T Coefficient type.
 * @param gNode Source node being processed (unused for counting).
 *
 * @note The thread index is obtained via mrcpp_get_thread_num().
 */
template <int D, typename T>
void OperatorStatistics::incrementGNodeCounters(const MWNode<D, T> &gNode) {
    int thread = mrcpp_get_thread_num();
    this->gCount[thread]++;
}

/**
 * @brief Increment the *f*-node application counters for the current thread.
 *
 * @tparam D Spatial dimension of the node.
 * @tparam T Coefficient type.
 * @param fNode Destination node to which an operator is applied.
 * @param ft    Destination component bitfield.
 * @param gt    Source component bitfield.
 *
 * @details
 * Increments:
 *  - Per-thread @c fCount.
 *  - Per-thread component histogram at entry (ft, gt).
 *  - Per-thread @c genCount if @c fNode.isGenNode() is true.
 */
template <int D, typename T>
void OperatorStatistics::incrementFNodeCounters(const MWNode<D, T> &fNode, int ft, int gt) {
    int thread = mrcpp_get_thread_num();
    this->fCount[thread]++;
    (*this->compCount[thread])(ft, gt) += 1;
    if (fNode.isGenNode()) { this->genCount[thread]++; }
}

/**
 * @brief Print a human-readable summary of accumulated totals.
 *
 * @param o Output stream.
 * @return Reference to @p o to allow chaining.
 *
 * @details
 * The output includes total counts for g-nodes, f-nodes, generalized nodes,
 * and the aggregated 8×8 (ft, gt) component histogram.
 */
std::ostream &OperatorStatistics::print(std::ostream &o) const {
    o << std::setw(8);
    o << "*OperatorFunc statistics: " << std::endl << std::endl;
    o << "  Total calculated gNodes      : " << this->totGCount << std::endl;
    o << "  Total applied fNodes         : " << this->totFCount << std::endl;
    o << "  Total applied genNodes       : " << this->totGenCount << std::endl << std::endl;
    o << "  By components:" << std::endl << *this->totCompCount << std::endl;
    return o;
}

/* ---- Explicit template instantiations for supported node types ---- */
template void OperatorStatistics::incrementFNodeCounters<1, double>(const MWNode<1, double> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<2, double>(const MWNode<2, double> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<3, double>(const MWNode<3, double> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<1, ComplexDouble>(const MWNode<1, ComplexDouble> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<2, ComplexDouble>(const MWNode<2, ComplexDouble> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<3, ComplexDouble>(const MWNode<3, ComplexDouble> &fNode, int ft, int gt);
template void OperatorStatistics::incrementGNodeCounters<1, double>(const MWNode<1, double> &gNode);
template void OperatorStatistics::incrementGNodeCounters<2, double>(const MWNode<2, double> &gNode);
template void OperatorStatistics::incrementGNodeCounters<3, double>(const MWNode<3, double> &gNode);
template void OperatorStatistics::incrementGNodeCounters<1, ComplexDouble>(const MWNode<1, ComplexDouble> &gNode);
template void OperatorStatistics::incrementGNodeCounters<2, ComplexDouble>(const MWNode<2, ComplexDouble> &gNode);
template void OperatorStatistics::incrementGNodeCounters<3, ComplexDouble>(const MWNode<3, ComplexDouble> &gNode);

} // namespace mrcpp
