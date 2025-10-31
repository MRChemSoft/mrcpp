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
 * @brief Adaptive multiwavelet convolution driver.
 *
 * @details
 * Declares @ref mrcpp::ConvolutionCalculator, a tree-walking calculator that
 * applies a (possibly non-local) @ref ConvolutionOperator to an input
 * @ref FunctionTree with adaptive precision control. The calculator
 * orchestrates band construction, per-node operator application, and optional
 * operator manipulation (e.g., unit-cell projections for periodic problems).
 */

#include "TreeCalculator.h"
#include "operators/OperatorStatistics.h"
#include "trees/FunctionTreeVector.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @class ConvolutionCalculator
 * @brief Performs adaptive convolution of a function tree with a convolution operator.
 *
 * @tparam D Spatial dimensionality (1–3).
 * @tparam T Coefficient scalar type (`double` or `ComplexDouble`).
 *
 * @details
 * The calculator traverses the output tree (owned by the base
 * @ref TreeCalculator) and, for each node, applies the convolution operator to
 * the relevant neighborhood (an operator *band*). Band sizes are derived from
 * the operator bandwidth and the current tree depth, and can be further tuned
 * by a user-supplied per-node precision function @ref setPrecFunction.
 *
 * The implementation records timing and operator statistics per band/component
 * to aid profiling, and can optionally manipulate the operator prior to
 * application (see @ref startManipulateOperator).
 *
 * ### Lifetime / ownership
 * - @ref ConvolutionCalculator does **not** own the operator nor the input
 *   function tree; it stores non-owning pointers.
 * - Timers and internal matrices are allocated and cleared by
 *   @ref initTimers / @ref clearTimers .
 */
template <int D, typename T>
class ConvolutionCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Construct a calculator for \f$ g = \mathcal{O}\{f\} \f$.
     *
     * @param p      Target accuracy (relative or absolute depending on usage).
     * @param o      Convolution operator to apply.
     * @param f      Input function tree \f$ f \f$.
     * @param depth  Maximum traversal depth for the output tree
     *               (defaults to @c MaxDepth for the MRA).
     *
     * @pre @p o and @p f must remain valid for the lifetime of the calculator.
     */
    ConvolutionCalculator(double p, ConvolutionOperator<D> &o, FunctionTree<D, T> &f, int depth = MaxDepth);

    /// @brief Destructor. Releases timers and internal band-size tables.
    ~ConvolutionCalculator() override;

    /**
     * @brief Produce the initial work vector of nodes for the output tree.
     *
     * @param tree Output tree that will receive the convolution result.
     * @return Pointer to a heap-allocated vector of nodes to start from.
     *
     * @details
     * The initial set typically includes end nodes (or generator nodes in
     * banded neighborhoods) where the operator action is non-zero.
     * The caller (base class) assumes ownership of the returned vector.
     */
    MWNodeVector<D, T>* getInitialWorkVector(MWTree<D, T> &tree) const override;

    /**
     * @brief Set a per-node precision function.
     *
     * @param prec_func A functor returning the local tolerance for a node index.
     *
     * @details
     * When provided, the calculator uses @p prec_func(idx) to refine the target
     * precision locally (e.g., tighter near singularities), typically in
     * conjunction with the global precision passed to the constructor.
     */
    void setPrecFunction(const std::function<double(const NodeIndex<D> &idx)> &prec_func) { this->precFunc = prec_func; }

    /**
     * @brief Enable operator manipulation prior to application.
     *
     * @param excUnit If `true`, manipulate on the unit cell (periodic contexts).
     *
     * @details
     * When enabled the operator may be preconditioned, symmetrized, or mapped
     * to a fundamental domain before application. Exact behavior depends on
     * the associated @ref ConvolutionOperator.
     */
    void startManipulateOperator(bool excUnit) {
        this->manipulateOperator = true;
        this->onUnitcell = excUnit;
    }

private:
    // ---- Configuration / inputs ------------------------------------------------

    /// @brief Maximum output depth to visit.
    int maxDepth;

    /// @brief Global target precision (interpreted by implementation).
    double prec;

    /// @brief Toggle for pre-application manipulation of the operator.
    bool manipulateOperator{false};

    /// @brief Toggle for unit-cell manipulation in periodic problems.
    bool onUnitcell{false};

    /// @brief Non-owning pointer to the convolution operator.
    ConvolutionOperator<D> *oper;

    /// @brief Non-owning pointer to the input function tree f(r).
    FunctionTree<D, T> *fTree;

    // ---- Instrumentation -------------------------------------------------------

    /// @brief Per-band timers for operator band building.
    std::vector<Timer*> band_t;

    /// @brief Per-band timers for the main convolution kernels.
    std::vector<Timer*> calc_t;

    /// @brief Per-band timers for norm/threshold checks.
    std::vector<Timer*> norm_t;

    /// @brief Aggregate operator statistics (bandwidths, touches, flops estimates).
    OperatorStatistics operStat;

    // ---- Band-size modeling ----------------------------------------------------

    /**
     * @brief Precomputed band-size factors per depth/component.
     *
     * @details
     * Each matrix has shape `(maxDepth+1) × nComp2`, where `nComp = 2^D`
     * and `nComp2 = nComp * nComp`. Linearized index
     * `k = gt * nComp + ft` maps from generator (`gt`) and father (`ft`)
     * component pairs to a band-size factor at a given depth.
     */
    std::vector<Eigen::MatrixXi*> bandSizes;

    /**
     * @brief Optional local precision override.
     *
     * @details
     * Defaults to a neutral functor returning 1.0. When set by
     * @ref setPrecFunction, it scales or replaces the global precision on a
     * per-node basis.
     */
    std::function<double(const NodeIndex<D> &idx)> precFunc = [](const NodeIndex<D> &idx) { return 1.0; };

    /// @brief Number of component blocks (2^D) in a multiwavelet tensor.
    static const int nComp = (1 << D);

    /// @brief Number of component-pair interactions ( (2^D) × (2^D) ).
    static const int nComp2 = (1 << D) * (1 << D);

    // ---- Band construction helpers --------------------------------------------

    /**
     * @brief Build an operator band (list of neighbor nodes) around @p gNode.
     *
     * @param gNode    Generator node (output-space anchor).
     * @param idx_band Output: collected node indices forming the band.
     * @return Heap-allocated node vector corresponding to @p idx_band.
     *
     * @details
     * The band is determined by the operator bandwidth at the scale of
     * @p gNode and the precomputed band-size factors. The returned vector
     * contains concrete node handles in traversal order.
     */
    MWNodeVector<D, T>* makeOperBand(const MWNode<D, T> &gNode, std::vector<NodeIndex<D>> &idx_band);

    /**
     * @brief Recursive fill of an operator band.
     *
     * @param band     Destination node vector to append to.
     * @param idx_band Indices to materialize.
     * @param idx      Current index under construction.
     * @param nbox     Periodic-box replication vector per dimension.
     * @param dim      Current dimension (0..D-1) being expanded.
     */
    void fillOperBand(MWNodeVector<D, T> *band, std::vector<NodeIndex<D>> &idx_band, NodeIndex<D> &idx, const int *nbox, int dim);

    // ---- Timing / statistics lifecycle -----------------------------------------

    /// @brief Allocate and start per-band timers.
    void initTimers();

    /// @brief Stop and free per-band timers.
    void clearTimers();

    /// @brief Print a compact timing breakdown per component/band.
    void printTimers() const;

    // ---- Band-size factors -----------------------------------------------------

    /// @brief Allocate @ref bandSizes tables.
    void initBandSizes();

    /**
     * @brief Lookup the band-size factor for a component pair at a depth.
     *
     * @param i      Which table (implementation-defined band decomposition).
     * @param depth  Tree depth.
     * @param os     Current operator-state (provides `gt` and `ft`).
     * @return Precomputed size factor.
     */
    int getBandSizeFactor(int i, int depth, const OperatorState<D, T> &os) const {
        int k = os.gt * this->nComp + os.ft;
        return (*this->bandSizes[i])(depth, k);
    }

    /**
     * @brief Compute the band-size factors for all component pairs at a depth.
     *
     * @param bs   Destination matrix (size `(maxDepth+1) × nComp2`).
     * @param depth Target depth to (re)compute.
     * @param bw   Operator bandwidth descriptor.
     */
    void calcBandSizeFactor(Eigen::MatrixXi &bs, int depth, const BandWidth &bw);

    // ---- Core calculation hooks (TreeCalculator overrides) ---------------------

    /**
     * @brief Compute the output contribution for a single node.
     *
     * @param node Output node to update.
     *
     * @details
     * Builds the relevant operator band around @p node, applies the operator to
     * the input tree restricted to that band, and accumulates the result into
     * @p node's coefficients. Precision is controlled by @ref prec and
     * @ref precFunc.
     */
    void calcNode(MWNode<D, T> &node) override;

    /**
     * @brief Post-processing after a full tree sweep.
     *
     * @details
     * Prints per-band timing information, clears timers, and re-initializes
     * the timing infrastructure for possible subsequent sweeps.
     */
    void postProcess() override {
        printTimers();
        clearTimers();
        initTimers();
    }

    // ---- Operator application kernels ------------------------------------------

    /**
     * @brief Apply a single operator component to the current band.
     *
     * @param os Operator state (component indices, buffers, thresholds, etc.).
     */
    void applyOperComp(OperatorState<D, T> &os);

    /**
     * @brief Apply the full operator (all components) for a given band index.
     *
     * @param i  Band table index / decomposition slot.
     * @param os Operator state for the current output node.
     */
    void applyOperator(int i, OperatorState<D, T> &os);

    /**
     * @brief Tensor-kernel variant of @ref applyOperComp (blocked/tensor form).
     *
     * @param os Operator state for the current output node.
     *
     * @details
     * May use batched/tensorized multiply-adds for better cache locality when
     * the component layout allows it.
     */
    void tensorApplyOperComp(OperatorState<D, T> &os);

    // ---- Tree maintenance -------------------------------------------------------

    /**
     * @brief Ensure parent nodes are materialized/touched before children writes.
     *
    * @param tree Output tree to touch/wake parents in.
     *
     * @details
     * Some backends require parent nodes to exist to safely commit child
     * contributions (e.g., for allocation, normalization, or boundary handling).
     */
    void touchParentNodes(MWTree<D, T> &tree) const;
};

} // namespace mrcpp