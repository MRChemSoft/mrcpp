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

#include "MWOperator.h"

namespace mrcpp {

/**
 * @class DerivativeOperator
 * @ingroup operators
 *
 * @brief Common base for derivative-type multiwavelet operators.
 *
 * @tparam D Spatial dimension of the operator (1, 2, or 3).
 *
 * @details
 * This abstract helper stores metadata and provides a thin interface for
 * operators that represent spatial derivatives in the multiwavelet (MW)
 * framework. It derives from @ref MWOperator and adds a single piece of
 * state: the derivative @ref order, which subclasses set appropriately
 * (e.g., 1 for first derivative, 2 for Laplacian-like second derivative
 * components, etc.).
 *
 * The constructor simply forwards the *scale window* to the base:
 * - @p root : the coarsest scale at which the operator is anchored,
 * - @p reach : the number of levels (half-width) the operator spans
 *   around @p root (default = 1).
 *
 * Concrete implementations such as @ref ABGVOperator and @ref BSOperator
 * specialize initialization, bandwidth, and stencil construction, while
 * reusing this small common interface.
 *
 * @see MWOperator, ABGVOperator, BSOperator
 */
template <int D> class DerivativeOperator : public MWOperator<D> {
public:
    /**
     * @brief Construct a derivative operator shell on a given scale window.
     *
     * @param mra   D-dimensional @ref MultiResolutionAnalysis that defines the
     *              domain and scaling basis for the operator.
     * @param root  Root scale (coarsest level) of the operator.
     * @param reach Scale reach around @p root (default = 1). A reach of @c r
     *              allows interaction across \f$ 2r+1 \f$ adjacent levels.
     *
     * @note This constructor does not assemble any stencil; subclasses call
     *       their own initialization routines and may update @ref order.
     */
    DerivativeOperator(const MultiResolutionAnalysis<D> &mra, int root, int reach = 1)
            : MWOperator<D>(mra, root, reach) {}

    DerivativeOperator(const DerivativeOperator &oper) = delete;             ///< Non-copyable
    DerivativeOperator &operator=(const DerivativeOperator &oper) = delete;  ///< Non-assignable
    ~DerivativeOperator() override = default;

    /**
     * @brief Return the derivative order encoded by this operator.
     * @returns Integer derivative order (1 by default; subclasses may set 2, 3, ...).
     */
    int getOrder() const { return order; }

protected:
    /** @brief Derivative order metadata (default = 1). Subclasses should set this. */
    int order{1};
};

} // namespace mrcpp
