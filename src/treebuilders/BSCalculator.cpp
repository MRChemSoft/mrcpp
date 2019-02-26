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

#include <fstream>

#include "MRCPP/config.h"

#include "BSCalculator.h"
#include "utils/Printer.h"

using Eigen::MatrixXd;

namespace mrcpp {

BSCalculator::BSCalculator(const ScalingBasis &basis, int n)
        : diff_order(n) {
    if (this->diff_order <= 0) NOT_IMPLEMENTED_ABORT;
    if (this->diff_order == 1) readSMatrix(basis, '1');
    if (this->diff_order == 2) readSMatrix(basis, '2');
    if (this->diff_order == 3) readSMatrix(basis, '3');
    if (this->diff_order >= 4) NOT_IMPLEMENTED_ABORT;
}

void BSCalculator::readSMatrix(const ScalingBasis &basis, char n) {
    std::string file;
    std::string path = MW_FILTER_DIR;
    if (basis.getScalingType() == Legendre) file = path + "/L_b-spline-deriv" + n + ".txt";
    if (basis.getScalingType() == Interpol) file = path + "/I_b-spline-deriv" + n + ".txt";
    if (basis.getScalingOrder() < 0) MSG_FATAL("Scaling order not supported");
    if (basis.getScalingOrder() > 20) MSG_FATAL("Scaling order not supported");

    std::ifstream ifs(file.c_str());
    if (not ifs) MSG_ERROR("Failed to open file: " << file);
    for (int kp1 = 2; kp1 < 21; kp1++) {
        std::string line;
        getline(ifs, line);
        std::istringstream iss(line);

        int order;
        iss >> order;
        if (order != kp1) MSG_FATAL("Orders no not match");

        MatrixXd data = MatrixXd::Zero(3 * kp1, kp1);
        for (int i = 0; i < 3 * kp1; i++) {
            getline(ifs, line);
            std::istringstream iss(line);
            for (int j = 0; j < kp1; j++) { iss >> data(i, j); }
        }
        if (kp1 == (basis.getScalingOrder() + 1)) {
            this->S_p1 = data.block(0 * kp1, 0, kp1, kp1);
            this->S_0 = data.block(1 * kp1, 0, kp1, kp1);
            this->S_m1 = data.block(2 * kp1, 0, kp1, kp1);
            break;
        }
    }
}

void BSCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();
    int np1 = node.getScale() + 1;
    int kp1 = node.getKp1();
    int kp1_d = node.getKp1_d();
    int l = node.getTranslation()[1] - node.getTranslation()[0];
    double two_np1 = std::pow(2.0, this->diff_order * np1);
    double *coefs = node.getCoefs();

    switch (l) {
        case 1:
            for (int i = 0; i < kp1; i++) {
                for (int j = 0; j < kp1; j++) {
                    int idx = i * kp1 + j;
                    coefs[1 * kp1_d + idx] = two_np1 * this->S_p1(i, j);
                }
            }
            break;
        case 0:
            for (int i = 0; i < kp1; i++) {
                for (int j = 0; j < kp1; j++) {
                    int idx = i * kp1 + j;
                    coefs[0 * kp1_d + idx] = two_np1 * this->S_0(i, j);
                    coefs[1 * kp1_d + idx] = two_np1 * this->S_m1(i, j);
                    coefs[2 * kp1_d + idx] = two_np1 * this->S_p1(i, j);
                    coefs[3 * kp1_d + idx] = two_np1 * this->S_0(i, j);
                }
            }
            break;
        case -1:
            for (int i = 0; i < kp1; i++) {
                for (int j = 0; j < kp1; j++) {
                    int idx = i * kp1 + j;
                    coefs[2 * kp1_d + idx] = two_np1 * this->S_m1(i, j);
                }
            }
            break;
        default:
            MSG_ERROR("This translation should not occour");
            break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

} // namespace mrcpp
