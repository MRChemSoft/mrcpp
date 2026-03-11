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
 * @brief Derivative calculator on multiresolution function trees.
 *
 * @details
 * Declares @ref mrcpp::DerivativeCalculator, a @ref TreeCalculator that applies a
 * directional differential operator to an input @ref FunctionTree and writes the
 * result into the calculator's target tree. The implementation constructs a
 * scale-aware operator “band” around each output node and evaluates the
 * derivative using tensorized component applications while optionally collecting
 * timing and bandwidth statistics.
 */

#include "TreeCalculator.h"
#include "operators/OperatorStatistics.h"

namespace mrcpp {

/**
 * @class DerivativeCalculator
 * @brief Applies a directional derivative operator to a function tree.
 *
 * @tparam D Spatial dimension of the tree (1–3 typical).
 * @tparam T Scalar coefficient type (e.g., `double`, `std::complex<double>`).
 *
 * @details
 * The calculator computes \f$ g = \partial_{x_{dir}}(f) \f$ where:
 *  - `dir` selects the Cartesian direction \f$0 \le dir < D\f$,
 *  - `oper` encapsulates the discretized derivative stencils/filters,
 *  - `fTree` is the source tree and the calculator’s target is the destination.
 *
 * The traversal is driven by the base @ref TreeCalculator; for each output
 * node, a localized operator band is constructed (see @ref makeOperBand) and
 * the operator is applied in a tensorized fashion (see
 * @ref tensorApplyOperComp). Optional timing is gathered per phase (band
 * building, application, norm updates) and summarized on completion.
 */
template <int D, typename T>
class DerivativeCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Construct a derivative calculator.
     *
     * @param dir Direction index of the derivative (0-based, `< D`).
     * @param o   Reference to the derivative operator to apply.
     * @param f   Reference to the **input/source** function tree.
     *
     * @note The destination/output tree is owned by the base
     *       @ref TreeCalculator (`this->outTree()`).
     */
    DerivativeCalculator(int dir, DerivativeOperator<D> &o, FunctionTree<D, T> &f);

    /// @brief Virtual destructor; prints and clears timers in @ref postProcess.
    ~DerivativeCalculator() override;

    /**
     * @brief Provide the initial work vector for traversal.
     *
     * @param tree Output tree on which work will be scheduled.
     * @return Pointer to a newly created vector of nodes to process first.
     *
     * @details
     * The default strategy is to populate the initial vector with those nodes
     * of @p tree that require operator application (implementation-dependent).
     */
    MWNodeVector<D, T> *getInitialWorkVector(MWTree<D, T> &tree) const override;

    /**
     * @brief Compute the derivative for a node pair.
     *
     * @param fNode Source node from @ref fTree (input).
     * @param gNode Destination node on the output tree (result of \f$\partial_{dir} f\f$).
     *
     * @details
     * Builds a local band around @p gNode, gathers contributions from @p fNode
     * within the operator bandwidth, then accumulates into @p gNode.
     */
    void calcNode(MWNode<D, T> &fNode, MWNode<D, T> &gNode);

private:
    // --- Configuration and inputs -------------------------------------------------

    int applyDir;                      ///< Direction index along which to differentiate.
    FunctionTree<D, T> *fTree{nullptr};///< Source function tree (input).
    DerivativeOperator<D> *oper{nullptr}; ///< Differential operator to apply.

    // --- Timing/statistics --------------------------------------------------------

    std::vector<Timer> band_t;  ///< Timers for band construction per depth or phase.
    std::vector<Timer> calc_t;  ///< Timers for operator application per depth or phase.
    std::vector<Timer> norm_t;  ///< Timers for norm/cleanup updates per depth or phase.
    OperatorStatistics operStat;///< Aggregate operator and bandwidth statistics.

    // --- Work preparation ---------------------------------------------------------

    /**
     * @brief Build the operator "band" (neighborhood) for an output node.
     *
     * @param gNode     Output (destination) node.
     * @param idx_band  Output: list of source node indices involved by bandwidth.
     * @return A vector of node pointers representing the band to be processed.
     *
     * @details
     * The band captures the set of input nodes that may contribute to @p gNode
     * under the derivative operator’s bandwidth model across scales and
     * spatial adjacency.
     */
    MWNodeVector<D, T> makeOperBand(const MWNode<D, T> &gNode, std::vector<NodeIndex<D>> &idx_band);

    /// @brief Initialize per-phase timers based on tree depth/layout.
    void initTimers();

    /// @brief Stop/clear all timers and release related resources.
    void clearTimers();

    /// @brief Print a concise timing summary and collected operator statistics.
    void printTimers() const;

    // --- TreeCalculator interface -------------------------------------------------

    /**
     * @brief Per-node callback from the traversal engine.
     *
     * @param node Output node to compute; pulls required input contributions.
     *
     * @details
     * For each output @p node, constructs the operator band (via
     * @ref makeOperBand), then delegates the actual stencil application to
     * @ref applyOperator / @ref tensorApplyOperComp. Norm and flag updates are
     * performed as needed.
     */
    void calcNode(MWNode<D, T> &node) override;

    /**
     * @brief Hook invoked after a traversal pass.
     *
     * @details Prints timing statistics, clears timers, and re-initializes
     * them to be ready for subsequent passes.
     */
    void postProcess() override {
        printTimers();
        clearTimers();
        initTimers();
    }

    // --- Operator application -----------------------------------------------------

    /**
     * @brief Apply the derivative operator to the current band/state.
     *
     * @param os Operator state for the current output node and component pair.
     *
     * @details
     * Chooses an application path depending on operator bandwidth and node
     * configuration, then accumulates results into the output node.
     */
    void applyOperator(OperatorState<D, T> &os);

    /**
     * @brief Specialized path for zero-bandwidth (local) derivative application.
     *
     * @param os Operator state for the current output node and component pair.
     *
     * @details
     * When the derivative is strictly local in the discretization (bandwidth 0),
     * this fast path avoids neighborhood assembly and directly applies the
     * local stencil.
     */
    void applyOperator_bw0(OperatorState<D, T> &os);

    /**
     * @brief Tensorized component application of the operator.
     *
     * @param os Operator state (contains gt/ft component ids, indices, buffers).
     *
     * @details
     * Performs dimension-wise application of the derivative operator using
     * separable tensor components, respecting the grid/scale layout in
     * @ref OperatorState.
     */
    void tensorApplyOperComp(OperatorState<D, T> &os);
};

} // namespace mrcpp