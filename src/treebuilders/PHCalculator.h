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
 * @brief 2D calculator that applies a three-point stencil in a scaling basis.
 *
 * @details
 * This component specializes the generic @ref TreeCalculator for 2D trees to
 * perform operations that can be expressed through **shifted overlap/stencil
 * matrices** of a scaling basis. The operator is represented by three
 * precomputed banded blocks
 * @f$S_{-1}, S_{0}, S_{+1}@f$, which correspond to interactions with the
 * left, center, and right neighbor positions along a given axis in the
 * multiresolution grid. These blocks are assembled from a provided
 * @ref ScalingBasis and then applied node-wise in @ref calcNode.
 *
 * Typical use-cases include discrete differential operators (e.g., first/second
 * derivatives) or narrow-band filters that can be written as a three-point
 * stencil in the scaling-function coefficient space. The effective stencil
 * "width" (derivative order) is indicated by @ref diff_order.
 */

#include <Eigen/Core>

#include "TreeCalculator.h"

namespace mrcpp {

/**
 * @class PHCalculator
 * @brief Applies a scaling-basis three-point stencil on a 2D multiresolution tree.
 *
 * @details
 * The calculator preloads three overlap/stencil matrices derived from a
 * @ref ScalingBasis:
 * - @ref S_m1 : shifted block for offset @f$-1@f$,
 * - @ref S_0  : central block for offset @f$0@f$,
 * - @ref S_p1 : shifted block for offset @f$+1@f$.
 *
 * During @ref calcNode, these blocks are combined to transform the node's
 * coefficients according to the chosen stencil (e.g., a centered finite
 * difference of order @ref diff_order). The exact algebra depends on the
 * basis; see the implementation of @ref readSMatrix.
 *
 * The class is marked @c final because it provides a complete node-level
 * implementation tailored to a 3-band stencil and is not intended for further
 * subclassing.
 */
class PHCalculator final : public TreeCalculator<2> {
public:
    /**
     * @brief Construct the calculator and preload stencil blocks.
     *
     * @param basis      Scaling basis from which the banded overlap/stencil
     *                   matrices are derived. The basis determines support,
     *                   moments, and thus the actual entries of the @f$S@f$
     *                   blocks.
     * @param n          Nominal stencil/derivative order (e.g., 1 for first
     *                   derivative, 2 for second). The value is stored in
     *                   @ref diff_order and may influence how the three
     *                   blocks are combined inside @ref calcNode.
     *
     * @post
     *  - @ref diff_order is set from @p n.
     *  - @ref S_m1, @ref S_0, @ref S_p1 are populated via @ref readSMatrix().
     */
    PHCalculator(const ScalingBasis &basis, int n);

private:
    /**
     * @brief Logical order of the stencil/differential operator to apply.
     *
     * @details
     * This does not change the *size* of the precomputed blocks, but can alter
     * how @ref S_m1, @ref S_0, and @ref S_p1 are linearly combined inside
     * @ref calcNode (e.g., centered first vs. second derivative weights).
     */
    const int diff_order;

    /// @name Precomputed banded overlap/stencil blocks
    /// @{
    /// Block corresponding to a shift of @f$-1@f$ grid unit(s).
    Eigen::MatrixXd S_m1;
    /// Central (unshifted) block.
    Eigen::MatrixXd S_0;
    /// Block corresponding to a shift of @f$+1@f$ grid unit(s).
    Eigen::MatrixXd S_p1;
    /// @}

    /**
     * @brief Node-level application of the three-point stencil.
     *
     * @param node Target 2D node whose coefficient vector is transformed
     *             in-place according to the assembled operator. The method
     *             is invoked by the traversal implemented in the base
     *             @ref TreeCalculator.
     *
     * @details
     * Conceptually, this computes (in scaling space)
     * @f[
     *   \mathbf{c}_{\text{out}} \;\leftarrow\;
     *   w_{-1}\,S_{-1}\,\mathbf{c}_{-1} \;+\;
     *   w_{0}\, S_{0}\, \mathbf{c}_{0} \;+\;
     *   w_{+1}\,S_{+1}\,\mathbf{c}_{+1},
     * @f]
     * where weights @f$w_{\cdot}@f$ depend on @ref diff_order and the chosen
     * discrete scheme, and @f$\mathbf{c}_{k}@f$ denotes the coefficient vector
     * at relative offset @f$k \in \{-1,0,+1\}@f$. The exact assembly is
     * implementation-specific and consistent with the supplied @ref ScalingBasis.
     */
    void calcNode(MWNode<2> &node);

    /**
     * @brief Populate one of the stencil matrices from the scaling basis.
     *
     * @param basis Scaling basis providing overlap and shift relations.
     * @param n     Selector for the matrix to load/build:
     *              chooses among @f$S_{-1}@f$, @f$S_{0}@f$, or @f$S_{+1}@f$.
     *              (The accepted values and encoding are implementation-defined,
     *              but conceptually map to offsets -1, 0, and +1.)
     *
     * @details
     * Extracts or assembles the band-limited matrix corresponding to a given
     * neighbor offset with respect to the scaling-function grid induced by
     * @p basis. The resulting block is stored into one of @ref S_m1, @ref S_0,
     * or @ref S_p1 depending on @p n.
     */
    void readSMatrix(const ScalingBasis &basis, char n);
};

} // namespace mrcpp