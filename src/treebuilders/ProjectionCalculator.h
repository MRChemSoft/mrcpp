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
 * @brief Node-wise projector from an analytic/representable function to a
 *        multiresolution function tree.
 *
 * @details
 * This calculator implements the core **projection kernel** that takes a
 * user-provided @ref RepresentableFunction and, node by node, produces
 * multiresolution coefficients for an output tree managed by the surrounding
 * @ref TreeCalculator pipeline.
 *
 * The calculator is agnostic to the scheduling/refinement policy: it only
 * defines how a *single* node is computed (see @ref calcNode). The driving
 * logic (initial work list, adaptors, termination) is handled by
 * @ref TreeCalculator and its collaborators.
 *
 * ### Coordinate scaling
 * A per-dimension @p scaling_factor is supplied at construction. It is applied
 * consistently during node evaluation to support anisotropic grid scalings,
 * unit conversions, or jacobian-like preconditioning. Use a vector of ones to
 * disable scaling.
 */

#include "TreeCalculator.h"

namespace mrcpp {

// Forward declaration; the concrete definition is provided by MRCPP headers.
template <int D, typename T> class RepresentableFunction;

/**
 * @class ProjectionCalculator
 * @brief Projects a @ref RepresentableFunction onto the active output tree.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type (e.g., `double`, `ComplexDouble`).
 *
 * @details
 * For each node requested by the @ref TreeCalculator scheduler, this class:
 *  1. Builds or fetches the corresponding input stencil/quadrature on the
 *     node’s support.
 *  2. Evaluates the supplied @ref RepresentableFunction at those points,
 *     applying the provided per-axis @ref scaling_factor.
 *  3. Computes the node’s scaling/wavelet coefficients, writes them to the
 *     output tree, and updates norms/metadata.
 *
 * The calculator itself does **not** decide where to refine; use an adaptor
 * (e.g., wavelet- or operator-based) with @ref TreeCalculator to drive
 * adaptivity from residuals or norm estimates.
 */
template <int D, typename T>
class ProjectionCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Construct a projector from an analytic/representable function.
     *
     * @param[in] inp_func        Function to be projected. The pointer is
     *                            stored and must remain valid for the lifetime
     *                            of the calculator.
     * @param[in] sf              Per-dimension scaling factors applied to local
     *                            coordinates before evaluating @p inp_func.
     *                            Set to `{1, …, 1}` for no scaling.
     *
     * @note The output target (tree) and traversal policy are provided by the
     *       owning @ref TreeCalculator context; this constructor only binds the
     *       callable and the evaluation scaling.
     */
    ProjectionCalculator(const RepresentableFunction<D, T> &inp_func,
                         const std::array<double, D> &sf)
            : func(&inp_func)
            , scaling_factor(sf) {}

private:
    /// Source function to be sampled on each node’s stencil.
    const RepresentableFunction<D, T> *func;

    /// Per-axis multiplicative coordinate scaling used during evaluation.
    const std::array<double, D> scaling_factor;

    /**
     * @brief Compute a single node of the output tree.
     *
     * @param[in,out] node Target node whose coefficients and norms are produced.
     *
     * @details
     * The typical implementation flow is:
     *  - derive the node’s physical coordinates from its @ref NodeIndex and
     *    apply @ref scaling_factor,
     *  - evaluate @ref func at the required sample points,
     *  - assemble scaling/wavelet coefficients and write them into @p node,
     *  - set “has coefficients” flags and update per-component norms.
     *
     * Thread-safety: the method only mutates @p node and uses read-only access
     * to @ref func and @ref scaling_factor, so it is safe under the usual
     * per-node parallel scheduling employed by @ref TreeCalculator.
     */
    void calcNode(MWNode<D, T> &node) override;
};

} // namespace mrcpp