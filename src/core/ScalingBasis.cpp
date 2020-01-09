/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "ScalingBasis.h"
#include "utils/Printer.h"

namespace mrcpp {

ScalingBasis::ScalingBasis(int k, int t)
        : type(t)
        , order(k) {
    if (this->order < 1) MSG_ABORT("Invalid scaling order");
    int q_order = getQuadratureOrder();
    this->quadVals = Eigen::MatrixXd::Zero(q_order, q_order);
    this->cvMap = Eigen::MatrixXd::Zero(q_order, q_order);
    this->vcMap = Eigen::MatrixXd::Zero(q_order, q_order);
}

void ScalingBasis::evalf(const double *r, Eigen::MatrixXd &vals) const {
    if (vals.rows() != this->funcs.size()) MSG_ERROR("Invalid argument");

    for (int d = 0; d < vals.cols(); d++) {
        for (int k = 0; k < vals.rows(); k++) { vals(k, d) = getFunc(k).evalf(r[d]); }
    }
}

const Eigen::MatrixXd &ScalingBasis::getCVMap(int operation) const {
    if (operation == Forward) {
        return this->cvMap;
    } else {
        return this->vcMap;
    }
}

bool ScalingBasis::operator==(const ScalingBasis &basis) const {
    if (this->type != basis.type) return false;
    if (this->order != basis.order) return false;
    return true;
}

bool ScalingBasis::operator!=(const ScalingBasis &basis) const {
    if (this->type != basis.type) return true;
    if (this->order != basis.order) return true;
    return false;
}

std::ostream &ScalingBasis::print(std::ostream &o) const {
    o << " polynomial order      : " << getScalingOrder() << std::endl;
    if (getScalingType() == Legendre) {
        o << " polynomial type       : Legendre";
    } else if (getScalingType() == Interpol) {
        o << " polynomial type       : Interpolating";
    } else {
        o << " polynomial type       : Unknown";
    }
    return o;
}

} // namespace mrcpp
