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

#include "IdentityConvolution.h"
#include "IdentityKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

/** @returns New IdentityConvolution object
 *  @param[in] mra: Which MRA the operator is defined
 *  @param[in] pr: Build precision, closeness to delta function
 *  @details This will project a kernel of a single gaussian with
 *  exponent sqrt(10/build_prec).
 */
template <int D>
IdentityConvolution<D>::IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double prec)
        : ConvolutionOperator<D>(mra) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 10.0;

    IdentityKernel<D> kernel(k_prec);
    this->initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

template <int D>
IdentityConvolution<D>::IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double prec, int root, int reach)
        : ConvolutionOperator<D>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 100.0;

    IdentityKernel<D> kernel(k_prec);
    this->initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

template class IdentityConvolution<1>;
template class IdentityConvolution<2>;
template class IdentityConvolution<3>;

} // namespace mrcpp