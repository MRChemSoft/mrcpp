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

#include "DerivativeOperator.h"

namespace mrcpp {

/** @class ABGVOperator
 *
 * @brief Derivative operator as defined by Alpert, Beylkin, Ginez and Vozovoi,
 * J Comp Phys 182, 149-190 (2002).
 *
 * NOTE: This is the recommended derivative operator for "cuspy" or discontinuous
 * functions. The BSOperator is recommended for smooth functions.
 */

template <int D> class ABGVOperator final : public DerivativeOperator<D> {
public:
    ABGVOperator(const MultiResolutionAnalysis<D> &mra, double a, double b);
    ABGVOperator(const ABGVOperator &oper) = delete;
    ABGVOperator &operator=(const ABGVOperator &oper) = delete;

protected:
    void initialize(double a, double b);
};

} // namespace mrcpp
