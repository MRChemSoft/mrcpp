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

/** @class PoissonOperator
 *
 * @brief Convolution with the Poisson Green's function kernel
 *
 * @details The Poisson kernel is approximated as a sum of Gaussian
 * functions in order to allow for separated application of the operator
 * in the Cartesian directions:
 * \f$ P(r-r') = \frac{1}{|r-r'|} \approx \sum_m^M \alpha_m e^{-\beta_m (r-r')^2} \f$
 */

class PoissonOperator final : public ConvolutionOperator<3> {
public:
    PoissonOperator(const MultiResolutionAnalysis<3> &mra, double prec);
    PoissonOperator(const MultiResolutionAnalysis<3> &mra, double prec, int root, int reach = 1);
    PoissonOperator(const PoissonOperator &oper) = delete;
    PoissonOperator &operator=(const PoissonOperator &oper) = delete;
};

} // namespace mrcpp
