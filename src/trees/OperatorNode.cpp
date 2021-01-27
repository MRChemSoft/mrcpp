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

#include "OperatorNode.h"
#include "NodeAllocator.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

void OperatorNode::dealloc() {
    int sIdx = this->serialIx;
    this->serialIx = -1;
    this->parentSerialIx = -1;
    this->childSerialIx = -1;
    this->tree->decrementNodeCount(this->getScale());
    this->tree->getNodeAllocator().deallocNodes(sIdx);
}

/** Calculate one specific component norm of the OperatorNode.
 *
 * OperatorNorms are defined as matrix 2-norms that are expensive to calculate.
 * Thus we calculate some cheaper upper bounds for this norm for thresholding.
 * First a simple vector norm, then a product of the 1- and infinity-norm. */
double OperatorNode::calcComponentNorm(int i) const {
    int depth = getDepth();
    double prec = getOperTree().getNormPrecision();
    double thrs = std::max(MachinePrec, prec / (8.0 * (1 << depth)));

    VectorXd coef_vec;
    this->getCoefs(coef_vec);

    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    const VectorXd &comp_vec = coef_vec.segment(i * kp1_d, kp1_d);
    const MatrixXd comp_mat = MatrixXd::Map(comp_vec.data(), kp1, kp1);

    double norm = 0.0;
    double vecNorm = comp_vec.norm();
    if (vecNorm > thrs) {
        double infNorm = math_utils::matrix_norm_inf(comp_mat);
        double oneNorm = math_utils::matrix_norm_1(comp_mat);
        if (std::sqrt(infNorm * oneNorm) > thrs) {
            double twoNorm = math_utils::matrix_norm_2(comp_mat);
            if (twoNorm > thrs) norm = twoNorm;
        }
    }
    return norm;
}

} // namespace mrcpp
