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

/**
 * @brief Construct a PHCalculator and load its stencil blocks.
 *
 * @param[in] basis      Scaling basis (type and order).
 * @param[in] n          Derivative order (1 or 2 supported).
 *
 * @throws NOT_IMPLEMENTED_ABORT if n <= 0 or n >= 3.
 *
 * @details
 * Based on \p n the constructor selects the corresponding set of PH derivative
 * blocks and loads them from disk. Supported files are:
 *  - Legendre:  L_ph_deriv_1.txt, L_ph_deriv_2.txt
 *  - Interpol:  I_ph_deriv_1.txt, I_ph_deriv_2.txt
 */
PHCalculator::PHCalculator(const ScalingBasis &basis, int n)
        : diff_order(n) {
    if (this->diff_order <= 0) NOT_IMPLEMENTED_ABORT;
    if (this->diff_order == 1) readSMatrix(basis, '1');
    if (this->diff_order == 2) readSMatrix(basis, '2');
    if (this->diff_order >= 3) NOT_IMPLEMENTED_ABORT;
}

/**
 * @brief Read PH derivative blocks from text files for the given basis.
 *
 * @param[in] basis  Scaling basis (provides type and order).
 * @param[in] n      Character '1' or '2' selecting derivative order.
 *
 * @details
 * The file format is:
 *  - First line per order k+1 (k+1 = 2..29): an integer "order" sentinel.
 *  - Followed by a 3*(k+1) by (k+1) table (row-major in the file) containing
 *    the vertically stacked blocks:
 *      [ S_{+1} ; S_{0} ; S_{-1} ]
 *
 * Only the block triple corresponding to the active basis order (kp1 = k+1)
 * is kept:
 *  - S_p1 = rows [0*kp1 .. 1*kp1-1]
 *  - S_0  = rows [1*kp1 .. 2*kp1-1]
 *  - S_m1 = rows [2*kp1 .. 3*kp1-1]
 *
 * @note
 *  - Supported scaling orders: 0..28 for Interpol/Legendre (kp1 in 2..29).
 *  - Files are discovered via details::find_filters().
 */
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

/**
 * @brief Fill 2D node coefficients by applying the PH derivative stencil.
 *
 * @param[in,out] node  2D MW node to populate (coefficients in scaling basis).
 *
 * @details
 * Let idx = (i0, i1) be the 2D node index and l = i1 - i0. Depending on l,
 * the appropriate neighbour coupling block is selected:
 *  - l = +1 : right neighbour uses S_p1
 *  - l =  0 : interior uses S_0 (diagonal) with off-diagonals S_{-1}, S_{+1}
 *  - l = -1 : left neighbour uses S_m1
 *
 * The coefficient tensor is laid out as 4 contiguous tiles for the four
 * tensor children, each of size kp1_d = (k+1)^2. For each tile we accumulate
 * the matrix-product contribution and rescale by
 *   two_np1 = 2^{diff_order * (scale+1)}.
 *
 * Finally, coefficients are transformed to MW (Compression), marked present,
 * and node norms are updated.
 *
 * @note
 *  - For periodic trees, indices outside the world box are ignored (no write).
 *  - The switch default does nothing by design (periodic handling is upstream).
 */
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
