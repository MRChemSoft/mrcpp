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

/**
 * @file PHCalculator.cpp
 * @brief Populate piecewise-homogeneous (PH) derivative stencil blocks for
 *        2D MW nodes and apply them as local operators.
 *
 * @details
 * The PH operator is applied on a 2D tensor-product node and uses three
 * nearest-neighbour coupling blocks along the refinement line:
 *  - S_m1 : block coupling to the left child   (l = -1)
 *  - S_0  : block coupling to the same child   (l =  0)
 *  - S_p1 : block coupling to the right child  (l = +1)
 *
 * For a node at scale j, the coefficients are scaled by 2^{diff_order*(j+1)}
 * to account for the dyadic scaling of derivatives in multiresolution analysis.
 *
 * The block matrices are read from precomputed text files (see @ref readSMatrix)
 * that depend on the scaling basis (Legendre or Interpolating) and the
 * derivative order (currently n = 1 or 2).
 */

#include "PHCalculator.h"

#include <fstream>

#include "MRCPP/config.h"

#include "utils/Printer.h"
#include "utils/details.h"

using Eigen::MatrixXd;

namespace mrcpp {

PHCalculator::PHCalculator(const ScalingBasis &basis, int n)
        : diff_order(n) {
    if (this->diff_order <= 0) NOT_IMPLEMENTED_ABORT;
    if (this->diff_order == 1) readSMatrix(basis, '1');
    if (this->diff_order == 2) readSMatrix(basis, '2');
    if (this->diff_order >= 3) NOT_IMPLEMENTED_ABORT;
}

void PHCalculator::readSMatrix(const ScalingBasis &basis, char n) {
    std::string file;
    std::string path = details::find_filters();

    if (basis.getScalingType() == Legendre) file = path + "/L_ph_deriv_" + n + ".txt";
    if (basis.getScalingType() == Interpol) file = path + "/I_ph_deriv_" + n + ".txt";
    if (basis.getScalingOrder() < 0) MSG_ABORT("Scaling order not supported");
    if (basis.getScalingOrder() > 29) MSG_ABORT("Scaling order not supported");

    std::ifstream ifs(file.c_str());
    if (not ifs) MSG_ERROR("Failed to open file: " << file);
    for (int kp1 = 2; kp1 < 30; kp1++) {
        std::string line;
        getline(ifs, line);
        std::istringstream iss(line);

        int order;
        iss >> order;
        if (order != kp1) MSG_ABORT("Orders no not match");

        MatrixXd data = MatrixXd::Zero(3 * kp1, kp1);
        for (int i = 0; i < 3 * kp1; i++) {
            getline(ifs, line);
            std::istringstream iss(line);
            for (int j = 0; j < kp1; j++) { iss >> data(i, j); }
        }
        if (kp1 == (basis.getScalingOrder() + 1)) {
            this->S_p1 = data.block(0 * kp1, 0, kp1, kp1);
            this->S_0  = data.block(1 * kp1, 0, kp1, kp1);
            this->S_m1 = data.block(2 * kp1, 0, kp1, kp1);
            break;
        }
    }
}

void PHCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();

    const auto &idx = node.getNodeIndex();
    int l = idx[1] - idx[0];          // neighbour offset along refinement line
    int np1 = idx.getScale() + 1;     // j+1, used in dyadic derivative scaling
    int kp1 = node.getKp1();          // k+1 (polynomial order + 1)
    int kp1_d = node.getKp1_d();      // (k+1)^2, tile size per child
    double two_np1 = std::pow(2.0, this->diff_order * np1);
    double *coefs = node.getCoefs();

    switch (l) {
        case 1: // right neighbour: only S_{+1} contributes
            for (int i = 0; i < kp1; i++) {
                for (int j = 0; j < kp1; j++) {
                    int idx = i * kp1 + j;
                    coefs[1 * kp1_d + idx] = two_np1 * this->S_p1(i, j);
                }
            }
            break;
        case 0: // interior: stencil spans S_0 (diagonal) and S_{-1}, S_{+1}
            for (int i = 0; i < kp1; i++) {
                for (int j = 0; j < kp1; j++) {
                    int idx = i * kp1 + j;
                    coefs[0 * kp1_d + idx] = two_np1 * this->S_0 (i, j);
                    coefs[1 * kp1_d + idx] = two_np1 * this->S_m1(i, j);
                    coefs[2 * kp1_d + idx] = two_np1 * this->S_p1(i, j);
                    coefs[3 * kp1_d + idx] = two_np1 * this->S_0 (i, j);
                }
            }
            break;
        case -1: // left neighbour: only S_{-1} contributes
            for (int i = 0; i < kp1; i++) {
                for (int j = 0; j < kp1; j++) {
                    int idx = i * kp1 + j;
                    coefs[2 * kp1_d + idx] = two_np1 * this->S_m1(i, j);
                }
            }
            break;
        default:
            // When periodic do nothing, else it should never end up here.
            break;
    }
    node.mwTransform(Compression); // convert to MW (wavelet) coefficients
    node.setHasCoefs();            // mark coefficients present
    node.calcNorms();              // update node/component norms
}

} // namespace mrcpp