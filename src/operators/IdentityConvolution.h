/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
#include "IdentityKernel.h"

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
    /** @returns New IdentityConvolution object
     *  @param[in] mra: Which MRA the operator is defined
     *  @param[in] pr: Build precision, closeness to delta function
     *  @details This will project a kernel of a single gaussian with
     *  exponent sqrt(10/build_prec).
     */
    IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double pr)
            : ConvolutionOperator<D>(mra, pr) {
        double epsilon = this->prec / 10.0;
        IdentityKernel identity_kernel(epsilon);
        this->initializeOperator(identity_kernel);
    }
    IdentityConvolution(const IdentityConvolution &oper) = delete;
    IdentityConvolution &operator=(const IdentityConvolution &oper) = delete;
};

} // namespace mrcpp
