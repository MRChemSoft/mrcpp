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
#include "DerivativeKernel.h"

namespace mrcpp {

/** @class DerivativeConvolution
 *
 * @brief Convolution with a derivative kernel
 *
 * @details The derivative kernel (derivative of Dirac's delta function) is
 * approximated by the derivative of a narrow Gaussian function:
 * \f$ D^x(r-r') = \frac{d}{dx}\delta(r-r') \approx \frac{d}{dx} \alpha e^{-\beta (r-r')^2} \f$
 */

template <int D> class DerivativeConvolution final : public ConvolutionOperator<D> {
public:
    /** @returns New DerivativeConvolution object
     *  @param[in] mra: Which MRA the operator is defined
     *  @param[in] pr: Build precision, closeness to delta function
     *  @details This will project a kernel of a single differentiated
     *  gaussian with exponent sqrt(10/build_prec).
     */
    DerivativeConvolution(const MultiResolutionAnalysis<D> &mra, double prec)
            : ConvolutionOperator<D>(mra, prec, prec / 10.0) {
        DerivativeKernel derivative_kernel(this->kern_prec);
        this->initializeOperator(derivative_kernel);
    }
    DerivativeConvolution(const MultiResolutionAnalysis<D> &mra, double prec, int root, int reach)
            : ConvolutionOperator<D>(mra, prec, prec / 100.0, root, reach) {
        DerivativeKernel derivative_kernel(this->kern_prec);
        this->initializeOperator(derivative_kernel);
    }
    DerivativeConvolution(const DerivativeConvolution &oper) = delete;
    DerivativeConvolution &operator=(const DerivativeConvolution &oper) = delete;
};

} // namespace mrcpp
