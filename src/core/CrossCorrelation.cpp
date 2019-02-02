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
 *
 *
 *  \date Jul 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include <fstream>

#include "MRCPP/config.h"
#include "MRCPP/constants.h"

#include "CrossCorrelation.h"
#include "utils/Printer.h"

using namespace Eigen;

namespace mrcpp {

std::string CrossCorrelation::default_ccc_lib = MW_FILTER_DIR;

CrossCorrelation::CrossCorrelation(int k, int t, const std::string &lib)
        : type(t)
        , order(k) {
    if (this->order < 1 or this->order > MaxOrder) { MSG_FATAL("Invalid cross correlation order: " << this->order); }
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    char *ep = getenv("MRCPP_FILTER_DIR");
    if (ep != nullptr) { default_ccc_lib = *ep; }
    int K = this->order + 1;
    setCCCPaths(lib);

    this->Left = MatrixXd(K * K, 2 * K);
    this->Right = MatrixXd(K * K, 2 * K);

    readCCCBin();
}

CrossCorrelation::CrossCorrelation(int t, const MatrixXd &L, const MatrixXd &R)
        : type(t)
        , order(L.cols() / 2 - 1) {
    if (this->order < 1 or this->order > MaxOrder) { MSG_FATAL("Invalid cross correlation order, " << this->order); }
    if (R.cols() != L.cols()) { MSG_FATAL("Right and Left cross correlation have different order!"); }
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    int K = this->order + 1;
    this->Left = MatrixXd(K * K, 2 * K);
    this->Right = MatrixXd(K * K, 2 * K);
    this->Left = L;
    this->Right = R;
}

void CrossCorrelation::setDefaultLibrary(const std::string &dir) {
    if (dir.empty()) { MSG_ERROR("No directory specified!"); }
    default_ccc_lib = dir;
}

void CrossCorrelation::setCCCPaths(const std::string &lib) {
    std::ostringstream oss;
    oss << this->order;
    std::string ordr = oss.str();
    std::string cclib;
    if (lib.empty()) {
        cclib = default_ccc_lib;
    } else {
        cclib = lib;
    }
    switch (this->type) {
        case (Interpol):
            this->L_path = cclib + "/I_c_left_" + ordr;
            this->R_path = cclib + "/I_c_right_" + ordr;
            break;
        case (Legendre):
            this->L_path = cclib + "/L_c_left_" + ordr;
            this->R_path = cclib + "/L_c_right_" + ordr;
            break;
        default:
            MSG_ERROR("Invalid CrossCorrelation type");
    }
}

void CrossCorrelation::readCCCBin() {
    std::ifstream L_fis(this->L_path.c_str(), std::ios::binary);
    std::ifstream R_fis(this->R_path.c_str(), std::ios::binary);

    if (not L_fis) { MSG_FATAL("Could not open cross correlation: " << this->L_path); }
    if (not R_fis) { MSG_FATAL("Could not open cross correlation: " << this->R_path); }

    int K = this->order + 1;
    double dL[2 * K];
    double dR[2 * K];
    for (int i = 0; i < K * K; i++) {
        L_fis.read((char *)dL, sizeof(double) * 2 * K);
        R_fis.read((char *)dR, sizeof(double) * 2 * K);
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
