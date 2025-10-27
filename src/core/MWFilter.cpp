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
 * Overview
 * --------
 * Implementation of the MWFilter class: a container for a 2K×2K multiwavelet
 * filter bank split into four K×K blocks (G0, G1, H0, H1) along with their
 * transposes. The filter bank supports two families (Interpol, Legendre) and
 * polynomial order 'order' (with K = order + 1).
 *
 * Block layout and semantics
 * --------------------------
 *   filter = [ G0  G1 ]   (top block-row: scaling/low-pass-like)
 *            [ H0  H1 ]   (bottom block-row: wavelet/high-pass-like)
 *
 * The precise interpretation (low-/high-pass) is family dependent, but the
 * layout is consistent across MRCPP. The class provides:
 *   - Loading G0 and H0 from binary files on disk.
 *   - Constructing G1 and H1 from symmetry relations (family-specific).
 *   - A full 2K×2K filter matrix 'filter' assembled from the four blocks.
 *   - Fast access to blocks and their transposes for compression/
 *     reconstruction phases of the multiresolution transform.
 *
 * File I/O conventions
 * --------------------
 *   - Files are discovered via details::find_filters() and named by family:
 *       Interpol: I_H0_<order>, I_G0_<order>
 *       Legendre: L_H0_<order>, L_G0_<order>
 *   - Format: raw binary doubles; K rows of K doubles each, row-major-by-row
 *     read in this implementation (one row per read).
 *   - Endianness and sizeof(double) must match the producing system.
 *
 * Symmetry completion
 * -------------------
 * Given H0 and G0 from disk, H1 and G1 are derived analytically:
 *   Interpol:
 *     G1(i,j) = (-1)^(i+K) * G0(i, K-j-1)
 *     H1(i,j) = H0(K-i-1, K-j-1)
 *   Legendre:
 *     G1(i,j) = (-1)^(i+j+K) * G0(i,j)
 *     H1(i,j) = (-1)^(i+j)   * H0(i,j)
 *
 * Transform directions
 * --------------------
 * - Reconstruction uses the blocks directly:      [H0 G0; H1 G1] in getSubFilter().
 * - Compression uses transposes of blocks:        [H0^T H1^T; G0^T G1^T].
 *   The mapping is encoded via getSubFilter(i, Compression/Reconstruction).
 *
 * Apply vs ApplyInverse
 * ---------------------
 * - apply(M/V):   multiplies by 'filter'          → reconstruction direction.
 * - applyInverse: multiplies by 'filter^T'        → compression direction.
 *   Both guard that the input vector/matrix has compatible row dimension 2K.
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

#include <fstream>
#include <iostream>

#include "MRCPP/config.h"
#include "MRCPP/constants.h"

#include "utils/Printer.h"
#include "utils/details.h"

using namespace Eigen;

namespace mrcpp {

/*
 * Constructor: MWFilter(int k, int t)
 * -----------------------------------
 * Build a filter bank of family 't' and order 'k'.
 * Steps:
 *   1) Validate order and type.
 *   2) Set file paths for H0/G0 based on family and order.
 *   3) Read H0 and G0 from disk; synthesize H1/G1 from symmetry rules.
 *   4) Assemble the full 2K×2K 'filter' matrix as [G0 G1; H0 H1].
 */
MWFilter::MWFilter(int k, int t)
        : type(t)
        , order(k) {
    if (this->order < 0 or this->order > MaxOrder) MSG_ABORT("Invalid filter order: " << this->order);
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    setFilterPaths(details::find_filters());

    generateBlocks();

    int K = this->order + 1;
    this->filter = MatrixXd::Zero(2 * K, 2 * K);
    this->filter << this->G0, this->G1, this->H0, this->H1;
}

/*
 * Constructor: MWFilter(int t, const MatrixXd& data)
 * --------------------------------------------------
 * Construct a filter bank directly from a provided 2K×2K matrix 'data'
 * (no disk I/O). The order is inferred as order = data.cols()/2 - 1.
 * After validation, the four K×K blocks and their transposes are extracted.
 */
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

/*
 * fillFilterBlocks()
 * ------------------
 * Slice the unified 2K×2K matrix 'filter' into the four K×K sub-blocks and
 * precompute their transposes. This is used after constructing from 'data'.
 */
void MWFilter::fillFilterBlocks() {
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

/*
 * getSubFilter(i, oper)
 * ---------------------
 * Retrieve one of the four K×K subfilters depending on transform 'oper':
 *   - Compression:    returns transposed blocks in order (H0^T, H1^T, G0^T, G1^T).
 *   - Reconstruction: returns direct blocks in order   (H0,   G0,   H1,   G1).
 * Index i ∈ {0,1,2,3} selects which block in the specified order.
 * Aborts on invalid index or oper.
 */
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

/*
 * Shorthand accessors for one direction only (avoid passing 'oper').
 * - getCompressionSubFilter(i): H0^T, H1^T, G0^T, G1^T (i=0..3)
 * - getReconstructionSubFilter(i): H0, G0, H1, G1      (i=0..3)
 */
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

/*
 * apply / applyInverse
 * --------------------
 * Multiply a vector/matrix by the filter or its transpose.
 * - apply(...)        : filter * data        (reconstruction direction)
 * - applyInverse(...) : filter^T * data      (compression direction)
 * Both validate row dimension matches the filter size (2K).
 */
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

/*
 * setFilterPaths(lib)
 * -------------------
 * Compose full file paths for H0 and G0 depending on family and order.
 * The prefix is 'I_' for Interpol and 'L_' for Legendre.
 */
void MWFilter::setFilterPaths(const std::string &lib) {
    switch (this->type) {
        case (Interpol):
            this->H_path = lib + "/I_H0_" + std::to_string(this->order);
            this->G_path = lib + "/I_G0_" + std::to_string(this->order);
            break;
        case (Legendre):
            this->H_path = lib + "/L_H0_" + std::to_string(this->order);
            this->G_path = lib + "/L_G0_" + std::to_string(this->order);
            break;
        default:
            MSG_ABORT("Invalid filter type " << this->type);
    }
}

/*
 * generateBlocks()
 * ----------------
 * Read H0 and G0 from binary files and synthesize H1/G1 from symmetry.
 * Finally, precompute all transposes.
 *
 * File format assumptions:
 *   - Each of H0 and G0 stores K rows; each row contains K doubles.
 *   - This function reads one row at a time into temporary buffers dH, dG.
 */
void MWFilter::generateBlocks() {
    std::ifstream H_fis(this->H_path.c_str(), std::ios::binary);
    std::ifstream G_fis(this->G_path.c_str(), std::ios::binary);

    if (H_fis.fail()) MSG_ABORT("Could not open filter: " << this->H_path);
    if (G_fis.fail()) MSG_ABORT("Could not open filter: " << this->G_path);

    int K = this->order + 1;

    double dH[K];
    double dG[K];
    /* read H0 and G0 from disk */
    this->G0 = Eigen::MatrixXd::Zero(K, K);
    this->H0 = Eigen::MatrixXd::Zero(K, K);
    for (int i = 0; i < K; i++) {
        H_fis.read((char *)dH, sizeof(double) * K);
        G_fis.read((char *)dG, sizeof(double) * K);
        for (int j = 0; j < K; j++) {
            this->G0(i, j) = dG[j]; // G0
            this->H0(i, j) = dH[j]; // H0
        }
    }
    G_fis.close();
    H_fis.close();

    /* fill H1 and G1 according to symmetry */
    this->G1 = Eigen::MatrixXd::Zero(K, K);
    this->H1 = Eigen::MatrixXd::Zero(K, K);
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