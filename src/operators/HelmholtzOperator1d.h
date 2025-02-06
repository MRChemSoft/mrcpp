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

/** @class HelmholtzOperator1d
 *
 * @brief Convolution with the Helmholtz Green's function kernel for the one-dimensional case
 *
 * @details The 1d Helmholtz kernel is a single Slater-type function
 * 
 * \f$ H(r-r') = \frac{1}{2k} e^{-k}{|x-x'|} \f$
 */


class HelmholtzOperator1d final : public ConvolutionOperator<3> {
public:
    HelmholtzOperator1d(const MultiResolutionAnalysis<3> &mra, double k, double prec);
    HelmholtzOperator1d(const MultiResolutionAnalysis<3> &mra, double k, double prec, int root, int reach = 1);
    HelmholtzOperator1d(const HelmholtzOperator1d &oper) = delete;
    HelmholtzOperator1d &operator=(const HelmholtzOperator1d &oper) = delete;
};

} // namespace mrcpp

