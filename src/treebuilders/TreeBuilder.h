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

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @class TreeBuilder
 * @brief Orchestrates adaptive construction and refinement of @ref MWTree objects.
 *
 * @tparam D Spatial dimension of the tree.
 * @tparam T Coefficient value type (e.g., `double`, `ComplexDouble`).
 *
 * @details
 * `TreeBuilder` coordinates three roles during adaptive computations:
 *  - a **calculator** (@ref TreeCalculator) that evaluates node data
 *    (coefficients, norms, metadata) on the current grid,
 *  - an **adaptor** (@ref TreeAdaptor) that decides which nodes should
 *    be refined (split),
 *  - the **tree** (@ref MWTree) that stores topology and coefficients.
 *
 * A typical adaptive build loop is:
 *  1. @ref calc to populate coefficients/norms on the current grid,
 *  2. @ref split to refine nodes selected by the adaptor,
 *  3. repeat (1–2) until no more splits occur or `maxIter` is reached.
 *
 * Some calculators maintain internal statistics/timers and may need a final
 * post-processing step; @ref build calls into the calculator appropriately.
 */
template <int D, typename T>
class TreeBuilder final {
public:
/**
 * @brief Adaptive build: iterate (calc → split) up to @p maxIter times.
 *
 * @param[in,out] tree       Target tree to (re)build/refine.
 * @param[in,out] calculator Calculator used to fill coefficients/norms on the current grid.
 * @param[in,out] adaptor    Refinement policy deciding which nodes to split.
 * @param[in]     maxIter    Upper bound on calc/split passes (use a small integer; non-positive means 0 passes).
 *
 * @details
 * The method performs:
 *  - an initial calc pass,
 *  - up to @p maxIter refinement passes, each performing:
 *      - split(tree, adaptor, passCoefs=true)
 *      - calc(tree, calculator)
 *  - any calculator post-processing hooks.
 *
 * Implementations typically stop early when split returns 0 (no new nodes).
 */
    void build(MWTree<D, T> &tree,
               TreeCalculator<D, T> &calculator,
               TreeAdaptor<D, T> &adaptor,
               int maxIter) const;

    /**
     * @brief Clear node data in @p tree using the provided @p calculator policy.
     *
     * @param[in,out] tree       Tree whose nodes should be cleared.
     * @param[in,out] calculator Calculator that defines how to reset per-node state.
     *
     * @details
     * Resets coefficient flags and cached norms to a consistent "empty" state.
     * This is useful before reusing a tree structure for another computation.
     */
    void clear(MWTree<D, T> &tree, TreeCalculator<D, T> &calculator) const;

    /**
     * @brief Compute/refresh coefficients and norms on the current grid.
     *
     * @param[in,out] tree       Tree to evaluate.
     * @param[in,out] calculator Calculator that implements per-node computation.
     *
     * @details
     * Traverses the active nodes (calculator-dependent strategy) and ensures
     * each leaf has consistent coefficients (scaling/wavelet) and derived norms.
     */
    void calc(MWTree<D, T> &tree, TreeCalculator<D, T> &calculator) const;

    /**
     * @brief Refine the tree topology according to @p adaptor policy.
     *
     * @param[in,out] tree    Tree subject to refinement.
     * @param[in,out] adaptor Adaptor deciding which nodes to split.
     * @param[in]     passCoefs
     *   If `true`, propagate or initialize child coefficients immediately
     *   (calculator-dependent behavior); if `false`, only topology changes
     *   are performed and coefficients are left for a subsequent @ref calc.
     *
     * @return Number of **new nodes** created (sum of all children inserted).
     *
     * @details
     * The method collects candidate leaves, applies the adaptor’s
     * `splitNodeVector`, and updates the tree topology. Implementations may
     * perform light-weight coefficient seeding when `passCoefs==true` to
     * improve the next calculation pass.
     */
    int split(MWTree<D, T> &tree, TreeAdaptor<D, T> &adaptor, bool passCoefs) const;

private:
    /**
     * @brief Aggregate the total scaling-norm over a set of nodes.
     *
     * @param vec Vector of node pointers to be reduced.
     * @return Sum (or calculator-defined aggregation) of scaling coefficients' norm.
     *
     * @details
     * Utility used by build loops for convergence checks and diagnostics.
     * The precise definition of “scaling norm” follows the node’s basis.
     */
    double calcScalingNorm(const MWNodeVector<D, T> &vec) const;

    /**
     * @brief Aggregate the total wavelet-norm over a set of nodes.
     *
     * @param vec Vector of node pointers to be reduced.
     * @return Sum (or calculator-defined aggregation) of wavelet coefficients' norm.
     *
     * @details
     * Utility used by build loops for refinement heuristics and stopping criteria.
     * The precise definition of “wavelet norm” follows the node’s basis.
     */
    double calcWaveletNorm(const MWNodeVector<D, T> &vec) const;
};

} // namespace mrcpp