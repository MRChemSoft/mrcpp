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
 * @file PHOperator.h
 * @brief Declaration of a Holoborodko-style smoothing derivative operator.
 *
 * @details
 * This header declares @ref mrcpp::PHOperator, a lightweight derivative operator
 * constructed from the smooth, low-noise differentiators introduced by
 * Pavel Holoborodko (see
 * <a href="http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/">
 * reference link</a>).
 *
 * The operator is assembled in the multiwavelet framework and is intended
 * primarily for experimentation/validation with smoothing differentiators.
 * For robust production work:
 * - use @ref mrcpp::ABGVOperator for functions with cusps/discontinuities,
 * - or @ref mrcpp::BSOperator for sufficiently smooth functions.
 *
 * @see mrcpp::DerivativeOperator, mrcpp::ABGVOperator, mrcpp::BSOperator
 */

#pragma once

#include "DerivativeOperator.h"

namespace mrcpp {

/**
 * @class PHOperator
 * @ingroup operators
 *
 * @brief Derivative operator based on Holoborodko’s smooth, low-noise differentiators.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 *
 * @details
 * This class derives from @ref DerivativeOperator and provides a separable,
 * single-component derivative approximation whose stencil is defined by the
 * Holoborodko differentiators. Internally, the concrete operator blocks are
 * produced by a PH-specific calculator and stored in an @ref OperatorTree.
 *
 * @note This is **not** the recommended operator for general calculations. Prefer
 * @ref ABGVOperator for non-smooth data and @ref BSOperator for smooth data.
 */
template <int D> class PHOperator final : public DerivativeOperator<D> {
public:
    /**
     * @brief Construct a PH-based derivative operator.
     *
     * @param mra   MultiResolutionAnalysis defining the domain and basis.
     * @param order Derivative order (typically 1 or 2).
     *
     * @warning Orders beyond those implemented by the underlying calculator
     *          are not supported.
     */
    PHOperator(const MultiResolutionAnalysis<D> &mra, int order);

    PHOperator(const PHOperator &oper) = delete;            ///< Non-copyable
    PHOperator &operator=(const PHOperator &oper) = delete; ///< Non-assignable

protected:
    /**
     * @brief Build and cache the internal operator representation.
     *
     * @details
     * Creates the PH calculator for the current scaling basis and requested order,
     * assembles an @ref OperatorTree with bandwidth control, transforms it to the
     * multiwavelet domain, and initializes the separable expansion.
     */
    void initialize();
};

} // namespace mrcpp