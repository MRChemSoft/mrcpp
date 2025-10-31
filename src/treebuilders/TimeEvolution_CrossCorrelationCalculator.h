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
 * @brief Time-evolution calculator based on cross-correlation kernels.
 *
 * @details
 * This header declares a node-local calculator that evaluates contributions
 * required for (imaginary or real time) Schrödinger evolution using
 * precomputed *J-power* integrals and a cross-correlation driver.
 *
 * Conceptually, for a two-dimensional function tree \f$f(\mathbf r)\f$ the
 * calculator applies (per node) a correlation-type update of the form
 * \f[
 *   g(\mathbf r) \;=\; \big(K * f\big)(\mathbf r)
 *   \;=\; \int_{\mathbb R^2} K(\mathbf r - \mathbf r')\, f(\mathbf r')\, d\mathbf r' ,
 * \f]
 * where the kernel \f$K\f$ and its (power) moments are provided through
 * the cross-correlation infrastructure and the \c JpowerIntegrals table.
 *
 * The boolean switch #imaginary selects which component of the complex-valued
 * kernel (or of the assembled integral) is used:
 * - `imaginary == false` → **real** part;
 * - `imaginary == true`  → **imaginary** part.
 *
 * The class is invoked by the tree execution engine (see @ref TreeCalculator)
 * and operates independently on each @ref MWNode.
 */

#include "TreeCalculator.h"
#include "core/CrossCorrelationCache.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"
#include "functions/JpowerIntegrals.h"

namespace mrcpp {

/**
 * @class TimeEvolution_CrossCorrelationCalculator
 * @brief Node calculator for Schrödinger time evolution via cross-correlation.
 *
 * @details
 * This calculator evaluates nodewise contributions needed for real- or
 * imaginary-time propagation in 2D using a cross-correlation representation
 * of the evolution operator. Precomputed integrals of the form
 * \f$ J_m = \int x^m\,K(x)\,dx \f$ (and higher-dimensional analogs) are
 * supplied through a map of @ref JpowerIntegrals instances indexed by
 * the power/order.
 *
 * ### Responsibilities
 * - Pull required cross-correlation data from a
 *   @ref SchrodingerEvolution_CrossCorrelation instance.
 * - Select **real** or **imaginary** contribution according to #imaginary.
 * - Assemble the per-node update and write the result to the output node.
 *
 * ### Threading / Parallelism
 * The class itself holds only non-owning pointers and simple references
 * to shared, read-only tables. It is thus re-entrant across nodes.
 * Synchronization and scheduling are handled at the @ref TreeCalculator layer.
 *
 * @note All pointer members are **non-owning**; the caller must ensure they
 *       remain valid for the lifetime of the calculator.
 */
class TimeEvolution_CrossCorrelationCalculator final : public TreeCalculator<2> {
public:
    /**
     * @brief Construct the calculator with auxiliary integral tables and a driver.
     *
     * @param[in] J
     *   Map from power/order (e.g., \f$m\f$) to the corresponding
     *   @ref JpowerIntegrals table. The calculator does **not** take ownership.
     * @param[in] cross_correlation
     *   Pointer to a @ref SchrodingerEvolution_CrossCorrelation driver that
     *   exposes kernel accessors / caches needed to assemble the correlation
     *   at node level. Non-owning.
     * @param[in] imaginary
     *   If `true`, use the **imaginary part** of the accumulated contribution;
     *   otherwise use the **real part**.
     *
     * @warning The map and the driver pointer must outlive this calculator.
     */
    TimeEvolution_CrossCorrelationCalculator(std::map<int, JpowerIntegrals *> &J,
                                             SchrodingerEvolution_CrossCorrelation *cross_correlation,
                                             bool imaginary)
            : J_power_inetgarls(J)
            , cross_correlation(cross_correlation)
            , imaginary(imaginary) {}

    /**
     * @brief Compute the contribution for one output node.
     *
     * @param[in,out] node
     *   The node to be written. The implementation typically:
     *   1) gathers the necessary kernel moments / cache entries,
     *   2) accumulates the cross-correlation at the node resolution,
     *   3) commits coefficients and refreshes norms/flags.
     *
     * @note The exact algebra (e.g., reconstruction/compression steps) is
     *       implemented in the corresponding source file.
     */
    void calcNode(MWNode<2> &node) override;

    /**
     * @brief Apply the cross-correlation operator at the granularity of a single node.
     *
     * @param[in,out] node The node to which the operator is applied.
     *
     * @details
     * This helper encapsulates the node-local application of the correlation
     * kernel using the caches provided by #cross_correlation and the moment
     * tables from #J_power_inetgarls. The #imaginary flag governs whether
     * the real or imaginary component of the final integral is extracted.
     *
     * @see calcNode
     */
    void applyCcc(MWNode<2> &node);

    // ---------------------------------------------------------------------
    // Public state (non-owning) — kept public to match existing interfaces.
    // ---------------------------------------------------------------------

    /**
     * @brief Precomputed kernel moment/integral tables, indexed by power.
     *
     * @note Non-owning pointers; the map must remain valid externally.
     */
    std::map<int, JpowerIntegrals *> J_power_inetgarls;

    /**
     * @brief Cross-correlation driver (non-owning).
     *
     * Provides access to kernel caches and auxiliary data needed to assemble
     * the correlation at a given node.
     */
    SchrodingerEvolution_CrossCorrelation *cross_correlation;

    /**
     * @brief Component selector for complex contributions.
     *
     * If `false`, the **real** part is used; if `true`, the **imaginary** part.
     */
    bool imaginary;
};

} // namespace mrcpp