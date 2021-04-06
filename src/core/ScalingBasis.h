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

#include <vector>

#include "MRCPP/constants.h"

#include "functions/Polynomial.h"

namespace mrcpp {

class ScalingBasis {
public:
    ScalingBasis(int k, int t);
    virtual ~ScalingBasis() = default;

    void evalf(const double *r, Eigen::MatrixXd &vals) const;

    Polynomial &getFunc(int k) { return this->funcs[k]; }
    const Polynomial &getFunc(int k) const { return this->funcs[k]; }

    int getScalingType() const { return this->type; }
    int getScalingOrder() const { return this->order; }
    int getQuadratureOrder() const { return this->order + 1; }

    const Eigen::MatrixXd &getQuadratureValues() const { return this->quadVals; }
    const Eigen::MatrixXd &getCVMap(int operation) const;

    bool operator==(const ScalingBasis &basis) const;
    bool operator!=(const ScalingBasis &basis) const;

    friend std::ostream &operator<<(std::ostream &o, const ScalingBasis &bas) { return bas.print(o); }

protected:
    const int type;
    const int order;
    Eigen::MatrixXd quadVals; // function values at quadrature pts
    Eigen::MatrixXd cvMap;    // coef-value transformation matrix
    Eigen::MatrixXd vcMap;    // value-coef transformation matrix
    std::vector<Polynomial> funcs;

    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp
