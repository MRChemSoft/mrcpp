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

#include "functions/GaussExp.h"
#include "functions/GaussFunc.h"

namespace mrcpp {

/** @class HeatKernel.
 *
 * @brief Heat kernel in \f$ \mathbb R^D \f$.
 *
 * @details In $\mathbb R^D$ the heat kernel has the form
 * \f[
 *  K_t(x)
 *  =
 *   \frac 1{ (4 \pi t)^{D/2} }
 *   \exp
 *   \left(
 *       - \frac{ |x|^2 }{4t}
 *   \right)
 *   , \quad
 *   x \in \mathbb R^D
 *   \text{ and }
 *   t > 0
 *   .
 * \f]
 * 
 */
template <int D> class HeatKernel final : public GaussExp<1> {
public:
    HeatKernel(double t)
            : GaussExp<1>() {
        double expo = 0.25 / t;
        double coef = std::pow(expo / mrcpp::pi, D / 2.0);
        GaussFunc<1> gFunc(expo, coef);
        this->append(gFunc);
    }
};

} // namespace mrcpp
