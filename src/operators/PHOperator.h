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

/** @class PHOperator
 *
 *  @brief Derivative operator based on the smoothing derivative of
 *  <a
 *  href="http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/">
 *  Pavel Holoborodko
 *  </a>.
 *
 *  NOTE: This is _not_ the recommended derivative operator for practial calculations, it's
 *  a proof-of-concept operator. Use the ABGVOperator for "cuspy" functions and the
 *  BSOperator for smooth functions.
 */

template <int D> class PHOperator final : public DerivativeOperator<D> {
public:
    PHOperator(const MultiResolutionAnalysis<D> &mra, int order);
    PHOperator(const PHOperator &oper) = delete;
    PHOperator &operator=(const PHOperator &oper) = delete;

protected:
    void initialize();
};

} // namespace mrcpp
