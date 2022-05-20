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

/** @class IdentityConvolution
 *
 * @brief Convolution with an identity kernel
 *
 * @details The identity kernel (Dirac's delta function) is approximated by a
 * narrow Gaussian function:
 * \f$ I(r-r') = \delta(r-r') \approx \alpha e^{-\beta (r-r')^2} \f$
 */

template <int D> class IdentityConvolution final : public ConvolutionOperator<D> {
public:
    IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double prec) : IdentityConvolution(mra, prec, mra.getRootScale(), -1) {}
    IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double prec, int root, int reach = -1);
    IdentityConvolution(const IdentityConvolution &oper) = delete;
    IdentityConvolution &operator=(const IdentityConvolution &oper) = delete;
};

} // namespace mrcpp
