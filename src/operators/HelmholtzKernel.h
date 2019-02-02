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

namespace mrcpp {

class HelmholtzKernel final : public GreensKernel {
public:
    HelmholtzKernel(double m, double eps, double r_min, double r_max)
            : GreensKernel(eps, r_min, r_max)
            , mu(m) {
        initializeKernel();
    }
    HelmholtzKernel(const HelmholtzKernel &kern) = delete;
    HelmholtzKernel &operator=(const HelmholtzKernel &kern) = delete;

protected:
    const double mu; /**< exponent */
    void initializeKernel() override;
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp
