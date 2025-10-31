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
 * @file BSCalculator.cpp
 * @brief Local stencil builder for smooth multiresolution derivative operators (“BS” family).
 *
 * @details
 * The **BSCalculator** assembles the *local* building blocks used by the smooth
 * derivative operator (see BSOperator). For a chosen scaling basis and derivative
 * order \( n\in\{1,2,3\} \), it loads three pretabulated coupling matrices
 * \f$S_{-1}, S_{0}, S_{+1}\f$ which represent the action of the derivative on a
 * 1D scaling block and its immediate neighbors (left, center, right) at a given scale.
 *
 * Source of the matrices:
 * - Files are looked up via `details::find_filters()`.
 * - Filenames depend on the scaling basis type and derivative order:
 *   - Legendre scaling: `L_b-spline-deriv{n}.txt`
 *   - Interpolating scaling: `I_b-spline-deriv{n}.txt`
 * - For each supported polynomial order `kp1 = 2..20`, the file stores a stacked
 *   3·kp1 × kp1 array that is split into the three kp1 × kp1 blocks
 *   \f$S_{+1}\f$, \f$S_{0}\f$, \f$S_{-1}\f$ (in that order).
 *
 * Application on a node:
 * - Given a 2D operator node (with index difference \f$\ell = i_1 - i_0 \in \{-1,0,+1\}\f$),
 *   BSCalculator writes the appropriate block(s) into the 2×2 corner layout of the node:
 *   - \f$\ell = -1\f$: left-neighbor coupling uses \f$S_{-1}\f$
 *   - \f$\ell = 0 \f$: center block uses \f$S_{0}\f$, off-diagonals use \f$S_{\pm 1}\f$
 *   - \f$\ell = +1\f$: right-neighbor coupling uses \f$S_{+1}\f$
 *
 * Scale factor:
 * - Derivatives scale as \f$2^{n\,(j+1)}\f$ where \f$n\f$ is the derivative order
 *   and \f$j+1\f$ is the node scale `np1`. The calculator multiplies all filled
 *   entries by \f$2^{n\,(j+1)}\f$.
 *
 * Limits and errors:
 * - Supported derivative orders: 1, 2, 3.
 * - Supported scaling orders: 1..20 (i.e., `kp1 = 2..21` in MRCPP terminology).
 * - On unsupported cases or missing files, the code aborts with a diagnostic.
 */

#include "BSCalculator.h"

#include <fstream>

#include "MRCPP/config.h"

#include "utils/Printer.h"
#include "utils/details.h"

using Eigen::MatrixXd;

namespace mrcpp {

/**
 * @brief Construct a BSCalculator and load derivative coupling blocks.
 *
 * @param basis Scaling basis (determines file family and polynomial order).
 * @param n     Derivative order (1, 2, or 3).
 *
 * @details
 * Dispatches to #readSMatrix to load \f$S_{-1}, S_{0}, S_{+1}\f$ for the given basis
 * and derivative order. Orders \f$n \ge 4\f$ are not implemented.
 *
 * @throws Aborts on unsupported derivative order, unsupported scaling order,
 *         or if the filter file cannot be opened.
 */
BSCalculator::BSCalculator(const ScalingBasis &basis, int n)
        : diff_order(n) {
    if (this->diff_order <= 0) NOT_IMPLEMENTED_ABORT;
    if (this->diff_order == 1) readSMatrix(basis, '1');
    if (this->diff_order == 2) readSMatrix(basis, '2');
    if (this->diff_order == 3) readSMatrix(basis, '3');
    if (this->diff_order >= 4) NOT_IMPLEMENTED_ABORT;
}

/**
 * @brief Load the pretabulated derivative coupling matrices from disk.
 *
 * @param basis Scaling basis (type and order).
 * @param n     Character identifying derivative order: '1', '2' or '3'.
 *
 * @details
 * - Chooses filename by basis type and derivative order.
 * - Iterates over the entries in the file for polynomial orders `kp1 = 2..20`
 *   until it matches the current basis order (`basis.getScalingOrder() + 1`).
 * - Splits the stacked 3·kp1 × kp1 array into three kp1 × kp1 blocks:
 *   \f$S_{+1}\f$, \f$S_{0}\f$, \f$S_{-1}\f$.
 *
 * File format expectations (per `kp1` section):
 * - First line: integer `order` (must equal `kp1`).
 * - Next 3·kp1 lines: kp1 numbers per line (row-major), forming the stacked matrix.
 *
 * @throws Aborts if:
 * - the file cannot be opened,
 * - the on-file order header does not match the expected `kp1`,
 * - the basis scaling order is unsupported.
 */
void BSCalculator::readSMatrix(const ScalingBasis &basis, char n) {
    std::string file;
    std::string path = details::find_filters();

    if (basis.getScalingType() == Legendre) file = path + "/L_b-spline-deriv" + n + ".txt";
    if (basis.getScalingType() == Interpol) file = path + "/I_b-spline-deriv" + n + ".txt";
    if (basis.getScalingOrder() < 0) MSG_ABORT("Scaling order not supported");
    if (basis.getScalingOrder() > 20) MSG_ABORT("Scaling order not supported");

    std::ifstream ifs(file.c_str());
    if (not ifs) MSG_ERROR("Failed to open file: " << file);
    for (int kp1 = 2; kp1 < 21; kp1++) {
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
 * @brief Populate a 2D operator node with the appropriate local derivative blocks.
 *
 * @param node Operator node to fill (corner layout, 2×2 logical structure).
 *
 * @details
 * Let \f$\ell = \text{idx}[1] - \text{idx}[0]\f$ denote the neighbor offset in the
 * second minus the first index direction. Depending on \f$\ell\f$, write the relevant
 * coupling block(s) into the node storage and multiply all entries by the scale factor
 * \f$2^{\,\text{diff\_order}\cdot (j+1)}\f$, where \f$j+1 = \text{idx.getScale()}+1\f$.
 *
 * Block placement (coefficient planes are enumerated in the code as 0,1,2,3):
 * - \f$\ell = +1\f$: only the “+1” plane is filled with \f$S_{+1}\f$.
 * - \f$\ell = 0 \f$: planes [0,1,2,3] are filled with \f$S_{0}, S_{-1}, S_{+1}, S_{0}\f$.
 * - \f$\ell = -1\f$: only the “+2” plane is filled with \f$S_{-1}\f$.
 *
 * After filling, the node is transformed (Compression), flagged as having coefficients,
 * and its norms are computed.
 */
void BSCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();

    const auto &idx = node.getNodeIndex();
    int l = idx[1] - idx[0];
    int np1 = idx.getScale() + 1;
    int kp1 = node.getKp1();
    int kp1_d = node.getKp1_d();
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
            // When periodic do nothing, else it should never end up here.
            break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

} // namespace mrcpp