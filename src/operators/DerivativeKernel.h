/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "GreensKernel.h"
#include "functions/GaussFunc.h"
#include "functions/GaussPoly.h"

namespace mrcpp {

class DerivativeKernel final : public GreensKernel {
public:
    DerivativeKernel(double eps)
            : GreensKernel(eps, -1.0, -1.0) {
        initializeKernel();
    }
    DerivativeKernel(const DerivativeKernel &kern) = delete;
    DerivativeKernel &operator=(const DerivativeKernel &kern) = delete;

protected:
    void initializeKernel() {
        double alpha = 1.0 / this->epsilon;
        double coef = std::pow(alpha / mrcpp::pi, 1.0 / 2.0);
        GaussFunc<1> g(alpha, coef);
        GaussPoly<1> dg = g.differentiate(0);
        this->append(dg);
    }

    std::ostream &print(std::ostream &o) const {
        o << " DerivativeKernel: " << std::endl;
        o << " epsilon:  " << this->epsilon << std::endl;
        return o;
    }
};

} // namespace mrcpp
