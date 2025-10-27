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
 *  File purpose (high level):
 *  --------------------------
 *  This implementation provides the CrossCorrelation class used to load and
 *  store *cross-correlation coefficient matrices* for multiwavelet filters.
 *  Two families of filters are supported (as encoded by `type`):
 *    - Interpolatory (prefix "I_")
 *    - Legendre    (prefix "L_")
 *
 *  Given an integer "order" k (poly order), we define K = k + 1. The class
 *  expects to find two binary files that contain the left and right cross
 *  correlation blocks:
 *      <lib>/<P>_c_left_<order>
 *      <lib>/<P>_c_right_<order>
 *  where <P> is "I" or "L" depending on the family. The directory <lib> is
 *  discovered via `details::find_filters()`.
 *
 *  Each file stores K*K rows, and each row contains 2*K doubles. The data are
 *  read into two Eigen matrices:
 *      Left  : (K*K) x (2K)
 *      Right : (K*K) x (2K)
 *
 *  Notes on indexing and sizes:
 *    - K = order + 1
 *    - The (K*K) rows represent a flattened 2D (i,j) index; i,j = 0..K-1.
 *    - Each row has 2K columns; the "2K" arises from the two-sided support
 *      of the correlation stencil (negative and positive offsets).
 *
 *  Error handling:
 *    - The code uses MRCPP's messaging macros (MSG_ABORT / MSG_ERROR) to
 *      report invalid input or missing files.
 *
 *  Endianness / portability:
 *    - Files are read as raw binary `double`. They must be produced on an
 *      architecture with compatible endianness and `double` layout.
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

// ----------------------------------------------------------------------------
// Constructor: CrossCorrelation(int k, int t)
// Creates an object for filter family `t` (see CrossCorrelation.h for the
// enum/type codes) and polynomial order `k`. It validates the order, validates
// the family, discovers the filter library directory, composes the filenames,
// and immediately loads the binary data into `Left` and `Right`.
// ----------------------------------------------------------------------------
CrossCorrelation::CrossCorrelation(int k, int t)
        : type(t)
        , order(k) {
    // Sanity check on order. `MaxOrder` is a library constant limiting k.
    if (this->order < 1 or this->order > MaxOrder) MSG_ABORT("Invalid cross correlation order: " << this->order);

    // Validate filter family (currently Interpol or Legendre are accepted).
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    // Locate the directory holding precomputed filter/correlation files.
    // `details::find_filters()` returns the absolute path to that directory.
    setCCCPaths(details::find_filters());

    // Load binary matrices Left and Right from disk into Eigen::MatrixXd.
    readCCCBin();
}

// ----------------------------------------------------------------------------
// Constructor: CrossCorrelation(int t, const MatrixXd& L, const MatrixXd& R)
// Directly construct a CrossCorrelation from matrices already in memory.
// `order` is inferred from the number of columns: 2K columns → K = order+1.
// Ensures the Left/Right shapes are compatible and the family type is valid.
// No file I/O is performed here.
// ----------------------------------------------------------------------------
CrossCorrelation::CrossCorrelation(int t, const MatrixXd &L, const MatrixXd &R)
        : type(t)
        , order(L.cols() / 2 - 1) {
    // Derive order from matrix width (2K columns → order = K - 1).
    if (this->order < 1 or this->order > MaxOrder) MSG_ABORT("Invalid cross correlation order, " << this->order);
    if (R.cols() != L.cols()) MSG_ABORT("Right and Left cross correlation have different order!");

    // Validate family.
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    // Shallow copies into class members (Eigen handles the allocation).
    this->Left = L;
    this->Right = R;
}

// ----------------------------------------------------------------------------
// setCCCPaths: Compose the on-disk file paths for the left/right matrices.
// Input: `lib` is the directory returned by `details::find_filters()`.
// The filenames follow the convention:
//   Interpol: I_c_left_<order>,  I_c_right_<order>
//   Legendre: L_c_left_<order>,  L_c_right_<order>
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// readCCCBin: Open the two binary files and load them into Eigen matrices.
// File structure:
//   - Let K = order + 1.
//   - Each file contains K*K consecutive rows.
//   - Each row stores 2*K doubles (contiguous), representing one stencil line.
// Post-processing:
//   - Any absolute value < MachinePrec is zeroed to improve sparsity/readability.
//   - Matrices are resized to (K*K) x (2*K).
// ----------------------------------------------------------------------------
void CrossCorrelation::readCCCBin() {
    // Open both files in binary mode; abort if either is missing.
    std::ifstream L_fis(this->L_path.c_str(), std::ios::binary);
    std::ifstream R_fis(this->R_path.c_str(), std::ios::binary);

    if (not L_fis) MSG_ABORT("Could not open cross correlation: " << this->L_path);
    if (not R_fis) MSG_ABORT("Could not open cross correlation: " << this->R_path);

    // Derive matrix dimensions from order.
    int K = this->order + 1;
    this->Left = MatrixXd::Zero(K * K, 2 * K);
    this->Right = MatrixXd::Zero(K * K, 2 * K);

    // Temporary row buffers for reading a single row (2K doubles) at a time.
    double dL[2 * K];
    double dR[2 * K];

    // Loop over all K*K rows and fill both Left and Right matrices.
    for (int i = 0; i < K * K; i++) {
        // Read one row for Left and one row for Right (raw binary doubles).
        L_fis.read((char *)dL, sizeof(double) * 2 * K);
        R_fis.read((char *)dR, sizeof(double) * 2 * K);

        // Copy into Eigen matrices with small-value cleanup.
        for (int j = 0; j < 2 * K; j++) {
            if (std::abs(dL[j]) < MachinePrec) dL[j] = 0.0; // numerical zeroing
            if (std::abs(dR[j]) < MachinePrec) dR[j] = 0.0;
            this->Left(i, j) = dL[j];
            this->Right(i, j) = dR[j];
        }
    }

    // Close streams (RAII would also close on destruction, but explicit is clear).
    L_fis.close();
    R_fis.close();
}

} // namespace mrcpp
