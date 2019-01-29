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

/*
 * \breif
 */

#pragma once

#include "functions/GaussExp.h"
#include "functions/Gaussian.h"

namespace mrcpp {

class GreensKernel : public GaussExp<1> {
public:
    GreensKernel(double eps, double r_min, double r_max)
            : GaussExp<1>()
            , epsilon(eps)
            , rMin(r_min)
            , rMax(r_max) {}
    GreensKernel(const GreensKernel &kern) = delete;
    GreensKernel &operator=(const GreensKernel &kern) = delete;
    ~GreensKernel() override = default;

    double getEpsilon() const { return this->epsilon; }
    double getRMin() const { return this->rMin; }
    double getRMax() const { return this->rMax; }

    void rescale(int d) {
        for (int i = 0; i < this->size(); i++) {
            Gaussian<1> &gauss = this->getFunc(i);
            double coef = std::pow(gauss.getCoef(), 1.0 / d);
            gauss.setCoef(coef);
        }
    }

    friend std::ostream &operator<<(std::ostream &o, const GreensKernel &kernel) { return kernel.print(o); }

protected:
    const double epsilon;
    const double rMin; /**< lower extreme */
    const double rMax; /**< upper extreme  (not used) */

    virtual void initializeKernel() = 0;

    virtual std::ostream &print(std::ostream &o) const = 0;
};

} // namespace mrcpp
