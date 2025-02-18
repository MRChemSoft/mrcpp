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

#include "HeatOperator.h"
#include "HeatKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

/** @brief Constructor  of the HeatOperator object
 *  @returns New HeatOperator object
 *  @param[in] mra: Which MRA the operator is defined
 *  @param[in] t: Time moment
 *  @param[in] prec: Build precision
 *  @details This will project a kernel of a single gaussian with
 *  exponent \f$ 1/(4t) \f$.
 *
 */
template <int D>
HeatOperator<D>::HeatOperator(const MultiResolutionAnalysis<D> &mra, double t, double prec)
        : ConvolutionOperator<D>(mra) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 10.0;

    HeatKernel<D> kernel(t);
    this->initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

/** @brief Constructor of the HeatOperator object in case of Periodic Boundary Conditions (PBC)
 *  @returns New IdentityConvolution object
 *  @param[in] mra: Which MRA the operator is defined
 *  @param[in] t: Time moment
 *  @param[in] prec: Build precision
 *  @param[in] root: root scale of operator.
 *  @param[in] reach: width at root scale (applies to periodic boundary conditions)
 *  @details This will project a kernel of a single gaussian with
 *           exponent \f$ 1/(4t) \f$.
 *           This version of the constructor
 *           is used for calculations within periodic boundary conditions (PBC).
 *           The \a root parameter is the coarsest negative scale at wich the operator
 *           is applied. The \a reach parameter is the bandwidth of the operator at
 *           the root scale. For details see \ref MWOperator
 *
 */
template <int D>
HeatOperator<D>::HeatOperator(const MultiResolutionAnalysis<D> &mra, double t, double prec, int root, int reach)
        : ConvolutionOperator<D>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 100.0;

    HeatKernel<D> kernel(t);
    this->initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

template class HeatOperator<1>;
template class HeatOperator<2>;
template class HeatOperator<3>;

} // namespace mrcpp
