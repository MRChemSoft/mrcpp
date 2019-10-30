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
 *  \date Jul 8, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "MWFilter.h"

#include "MRCPP/config.h"
#include "MRCPP/constants.h"

#include "coefficients/FilterData.h"
#include "utils/Printer.h"

using namespace Eigen;

namespace mrcpp {

MWFilter::MWFilter(int k, int t)
        : type(t)
        , order(k) {
    if (this->order < 1 or this->order > MaxOrder) MSG_ABORT("Invalid filter order: " << this->order);
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    generateBlocks();

    int K = this->order + 1;
    this->filter = MatrixXd::Zero(2 * K, 2 * K);
    this->filter << this->G0, this->G1, this->H0, this->H1;
}

MWFilter::MWFilter(int t, const MatrixXd &data) {
    this->type = t;
    this->order = data.cols() / 2 - 1;
    if (this->order < 0 or this->order > MaxOrder) MSG_ABORT("Invalid filter order " << this->order);
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << type);
    }

    this->filter = data;
    fillFilterBlocks();
}

void MWFilter::fillFilterBlocks() noexcept {
    int K = this->order + 1;
    this->G0 = this->filter.block(0, 0, K, K);
    this->G1 = this->filter.block(0, K, K, K);
    this->H0 = this->filter.block(K, 0, K, K);
    this->H1 = this->filter.block(K, K, K, K);
    this->G0t = this->G0.transpose();
    this->G1t = this->G1.transpose();
    this->H0t = this->H0.transpose();
    this->H1t = this->H1.transpose();
}

const MatrixXd &MWFilter::getSubFilter(int i, int oper) const {
    switch (oper) {
        case (Compression):
            switch (i) {
                case (0):
                    return this->H0t;
                case (1):
                    return this->H1t;
                case (2):
                    return this->G0t;
                case (3):
                    return this->G1t;
                default:
                    MSG_ABORT("Filter index out of bounds");
            }
            break;
        case (Reconstruction):
            switch (i) {
                case (0):
                    return this->H0;
                case (1):
                    return this->G0;
                case (2):
                    return this->H1;
                case (3):
                    return this->G1;
                default:
                    MSG_ABORT("Filter index out of bounds");
            }
            break;
        default:
            MSG_ABORT("Invalid wavelet transformation");
    }
}

const MatrixXd &MWFilter::getCompressionSubFilter(int i) const {
    switch (i) {
        case (0):
            return this->H0t;
        case (1):
            return this->H1t;
        case (2):
            return this->G0t;
        case (3):
            return this->G1t;
        default:
            MSG_ABORT("Filter index out of bounds");
    }
}

const MatrixXd &MWFilter::getReconstructionSubFilter(int i) const {
    switch (i) {
        case (0):
            return this->H0;
        case (1):
            return this->G0;
        case (2):
            return this->H1;
        case (3):
            return this->G1;
        default:
            MSG_ABORT("Filter index out of bounds");
    }
}

void MWFilter::apply(MatrixXd &data) const {
    if (data.rows() != this->filter.cols()) { INVALID_ARG_ABORT }
    data = this->filter * data;
}

void MWFilter::applyInverse(MatrixXd &data) const {
    if (data.rows() != this->filter.cols()) { INVALID_ARG_ABORT }
    data = this->filter.transpose() * data;
}

void MWFilter::apply(VectorXd &data) const {
    if (data.rows() != this->filter.cols()) { INVALID_ARG_ABORT }
    data = this->filter * data;
}

void MWFilter::applyInverse(VectorXd &data) const {
    if (data.rows() != this->filter.cols()) { INVALID_ARG_ABORT }
    data = this->filter.transpose() * data;
}

void MWFilter::generateBlocks() noexcept {
    this->G0 = get_G0(this->type, this->order);
    this->H0 = get_H0(this->type, this->order);

    int K = this->order + 1;
    this->G1 = Eigen::MatrixXd::Zero(K, K);
    this->H1 = Eigen::MatrixXd::Zero(K, K);

    /* fill H1 and G1 according to symmetry */
    switch (this->type) {
        case Interpol:
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    this->G1(i, j) = std::pow(-1.0, i + K) * this->G0(i, K - j - 1);
                    this->H1(i, j) = this->H0(K - i - 1, K - j - 1);
                }
            }
            break;
        case Legendre:
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    this->G1(i, j) = std::pow(-1.0, i + j + K) * this->G0(i, j);
                    this->H1(i, j) = std::pow(-1.0, i + j) * this->H0(i, j);
                }
            }
            break;
    }

    this->G0t = this->G0.transpose();
    this->G1t = this->G1.transpose();
    this->H0t = this->H0.transpose();
    this->H1t = this->H1.transpose();
}
} // namespace mrcpp
