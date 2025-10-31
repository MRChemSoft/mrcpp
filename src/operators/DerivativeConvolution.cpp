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

#include "DerivativeConvolution.h"
#include "DerivativeKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

template <int D>
DerivativeConvolution<D>::DerivativeConvolution(const MultiResolutionAnalysis<D> &mra, double prec)
        : ConvolutionOperator<D>(mra) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 10.0;

    DerivativeKernel<D> kernel(k_prec);
    this->initialize(kernel, k_prec, o_prec);

    Printer::setPrintLevel(oldlevel);
}

template <int D>
DerivativeConvolution<D>::DerivativeConvolution(const MultiResolutionAnalysis<D> &mra, double prec, int root, int reach)
        : ConvolutionOperator<D>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 100.0;

    DerivativeKernel<D> kernel(k_prec);
    this->initialize(kernel, k_prec, o_prec);

    Printer::setPrintLevel(oldlevel);
}

template class DerivativeConvolution<1>;
template class DerivativeConvolution<2>;
template class DerivativeConvolution<3>;

} // namespace mrcpp