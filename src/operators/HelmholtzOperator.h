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

/** @class HelmholtzOperator
 *
 * @brief Convolution with the Helmholtz Green's function kernel
 *
 * @details The Helmholtz kernel is approximated as a sum of gaussian functions
 * in order to allow for separated application of the operator in the Cartesian
 * directions:
 * \f$ H(r-r') = \frac{e^{-\mu|r-r'|}}{|r-r'|} \approx \sum_m^M \alpha_m e^{-\beta_m (r-r')^2} \f$
 */

class HelmholtzOperator final : public ConvolutionOperator<3> {
public:
    HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double m, double prec);
    HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double m, double prec, int root, int reach = 1);
    HelmholtzOperator(const HelmholtzOperator &oper) = delete;
    HelmholtzOperator &operator=(const HelmholtzOperator &oper) = delete;
};

} // namespace mrcpp
