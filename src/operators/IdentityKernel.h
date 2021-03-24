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

#include "GreensKernel.h"
#include "functions/GaussFunc.h"

namespace mrcpp {

class IdentityKernel final : public GreensKernel {
public:
    IdentityKernel(double eps)
            : GreensKernel(eps, -1.0, -1.0) {
        initializeKernel();
    }
    IdentityKernel(const IdentityKernel &kern) = delete;
    IdentityKernel &operator=(const IdentityKernel &kern) = delete;

protected:
    void initializeKernel() {
        double alpha = std::sqrt(1.0 / this->epsilon);
        double coef = std::pow(alpha / mrcpp::pi, 1.0 / 2.0);
        GaussFunc<1> gFunc(alpha, coef);
        this->append(gFunc);
    }
    std::ostream &print(std::ostream &o) const {
        o << "Kernel: " << std::endl;
        o << "epsilon:  " << this->epsilon << std::endl;
        o << "rMin:     " << this->rMin << std::endl;
        o << "rMax:     " << this->rMax << std::endl;
        return o;
    }
};

} // namespace mrcpp
