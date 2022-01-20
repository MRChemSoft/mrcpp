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
#include "functions/GaussPoly.h"

namespace mrcpp {

template <int D> class DerivativeKernel final : public GaussExp<1> {
public:
    DerivativeKernel(double epsilon)
            : GaussExp<1>() {
        double alpha = 1.0 / epsilon;
        double coef = std::pow(alpha / mrcpp::pi, D / 2.0);
        GaussFunc<1> g(alpha, coef);
        GaussPoly<1> dg = g.differentiate(0);
        this->append(dg);
    }
};

} // namespace mrcpp
