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

/*
 *
 *
 *  \date Jul 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "CrossCorrelation.h"

#include <fstream>

#include "MRCPP/config.h"
#include "MRCPP/constants.h"

#include "utils/Printer.h"
#include "utils/details.h"

using namespace Eigen;

namespace mrcpp {

CrossCorrelation::CrossCorrelation(int k, int t)
        : type(t)
        , order(k) {
    if (this->order < 1 or this->order > MaxOrder) MSG_ABORT("Invalid cross correlation order: " << this->order);
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    setCCCPaths(details::find_filters());

    readCCCBin();
}

CrossCorrelation::CrossCorrelation(int t, const MatrixXd &L, const MatrixXd &R)
        : type(t)
        , order(L.cols() / 2 - 1) {
    if (this->order < 1 or this->order > MaxOrder) MSG_ABORT("Invalid cross correlation order, " << this->order);
    if (R.cols() != L.cols()) MSG_ABORT("Right and Left cross correlation have different order!");
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    this->Left = L;
    this->Right = R;
}

void CrossCorrelation::setCCCPaths(const std::string &lib) {
    switch (this->type) {
        case (Interpol):
            this->L_path = lib + "/I_c_left_" + std::to_string(this->order);
            this->R_path = lib + "/I_c_right_" + std::to_string(this->order);
            break;
        case (Legendre):
            this->L_path = lib + "/L_c_left_" + std::to_string(this->order);
            this->R_path = lib + "/L_c_right_" + std::to_string(this->order);
            break;
        default:
            MSG_ERROR("Invalid CrossCorrelation type");
    }
}

void CrossCorrelation::readCCCBin() {
    std::ifstream L_fis(this->L_path.c_str(), std::ios::binary);
    std::ifstream R_fis(this->R_path.c_str(), std::ios::binary);

    if (not L_fis) MSG_ABORT("Could not open cross correlation: " << this->L_path);
    if (not R_fis) MSG_ABORT("Could not open cross correlation: " << this->R_path);

    int K = this->order + 1;
    this->Left = MatrixXd::Zero(K * K, 2 * K);
    this->Right = MatrixXd::Zero(K * K, 2 * K);
    std::vector<double> dL(2 * K);
    std::vector<double> dR(2 * K);
    for (int i = 0; i < K * K; i++) {
        L_fis.read(reinterpret_cast<char *>(dL.data()), sizeof(double) * 2 * K);
        R_fis.read(reinterpret_cast<char *>(dR.data()), sizeof(double) * 2 * K);
        for (int j = 0; j < 2 * K; j++) {
            if (std::abs(dL[j]) < MachinePrec) dL[j] = 0.0;
            if (std::abs(dR[j]) < MachinePrec) dR[j] = 0.0;
            this->Left(i, j) = dL[j];
            this->Right(i, j) = dR[j];
        }
    }

    L_fis.close();
    R_fis.close();
}

} // namespace mrcpp
