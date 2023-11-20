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

#include "ConvolutionOperator.h"

namespace mrcpp {

/** @class HeatOperator semigroup
 *
 * @brief Convolution with a heat kernel
 *
 * @details The exponential heat operator
 * \f$
 *   \exp \left( t \partial_x^2 \right)
 * \f$
 * can be regarded as a convolution operator in $L^2(\mathbb R)$
 * of the form
 * \[
 *   \exp \left( t \partial_x^2 \right)
 *   f(x)
 *   =
 *   \frac 1{ \sqrt{4 \pi t} }
 *   \int_\R
 *   \exp
 *   \left(
 *       - \frac{ (x - y)^2 }{4t}
 *   \right)
 *   f(y) dy
 *   , \quad
 *   t > 0
 *   .
 * \]
 * 
 */
template <int D> class HeatOperator final : public ConvolutionOperator<D> {
public:
    HeatOperator(const MultiResolutionAnalysis<D> &mra, double t, double prec);
    HeatOperator(const MultiResolutionAnalysis<D> &mra, double t, double prec, int root, int reach = 1);
    HeatOperator(const HeatOperator &oper) = delete;
    HeatOperator &operator=(const HeatOperator &oper) = delete;
};

} // namespace mrcpp
