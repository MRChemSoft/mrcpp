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

#include "HelmholtzOperator.h"
#include "HelmholtzKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

/** @returns New HelmholtzOperator object
 *  @param[in] mra: Which MRA the operator is defined
 *  @param[in] m: Exponential parameter of the operator
 *  @param[in] pr: Build precision, closeness to exp(-mu*r)/r
 *  @details This will construct a gaussian expansion to approximate
 *  exp(-mu*r)/r, and project each term into a one-dimensional MW operator.
 *  Subsequent application of this operator will apply each of the terms to
 *  the input function in all Cartesian directions.
 */
HelmholtzOperator::HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double mu, double prec)
        : ConvolutionOperator<3>(mra) {
    int oldlevel = Printer::setPrintLevel(0);

    double o_prec = prec;
    double k_prec = prec / 10.0;
    double r_min = this->MRA.calcMinDistance(k_prec);
    double r_max = this->MRA.calcMaxDistance();

    HelmholtzKernel kernel(mu, k_prec, r_min, r_max);
    initialize(kernel, k_prec, o_prec);

    Printer::setPrintLevel(oldlevel);
}

HelmholtzOperator::HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double mu, double prec, int root, int reach)
        : ConvolutionOperator<3>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    double o_prec = prec;
    double k_prec = prec / 100.0;
    double r_min = this->MRA.calcMinDistance(k_prec);
    double r_max = this->MRA.calcMaxDistance();

    // Adjust r_max for periodic world
    if (mra.getWorldBox().isPeriodic()) {
        auto rel_root = this->oper_root - this->MRA.getRootScale();
        r_max *= std::pow(2.0, -rel_root);
        r_max *= (2.0 * this->oper_reach) + 1.0;
    }
    HelmholtzKernel kernel(mu, k_prec, r_min, r_max);
    initialize(kernel, k_prec, o_prec);

    Printer::setPrintLevel(oldlevel);
}

} // namespace mrcpp
