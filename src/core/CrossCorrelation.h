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

#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

#include "MRCPP/config.h"

namespace mrcpp {

/**
 * @class CrossCorrelation
 * @brief Container/loader for multiwavelet cross-correlation coefficient tables.
 *
 * This class encapsulates the left/right cross-correlation matrices associated
 * with a chosen multiwavelet filter family and polynomial order.
 *
 *  • The filter "family" is identified by an integer @c type
 *    (e.g., Interpolatory or Legendre; concrete codes are defined elsewhere
 *    and validated in the implementation).
 *
 *  • The polynomial @c order is k ≥ 1. We use K = k + 1 for dimensions.
 *
 *  • Two dense matrices are held:
 *        Left  ∈ ℝ^{(K·K) × (2K)},
 *        Right ∈ ℝ^{(K·K) × (2K)}.
 *    Each row corresponds to a flattened (i,j) pair with i,j∈{0..K−1};
 *    each row stores a 2K-wide correlation stencil.
 *
 * Objects can be constructed by loading the binary coefficient files from disk
 * (constructor #1) or by adopting matrices already residing in memory
 * (constructor #2). Accessors expose the type/order and const references to
 * the matrices; there are no mutating public methods by design.
 *
 * Invariants (enforced in the implementation):
 *  - 1 ≤ order ≤ MaxOrder
 *  - Left.cols() == Right.cols() == 2K where K = order + 1
 *  - Left.rows() == Right.rows() == K*K
 *
 * Thread-safety: the class is a simple value holder once constructed.
 * Concurrent reads are safe; concurrent writes are not supported.
 */
class CrossCorrelation final {
public:
    /**
     * @brief Construct by loading coefficient tables from the filter library.
     *
     * The library path is discovered internally (see details::find_filters()).
     * Files are chosen based on @p type and @p k and read into #Left/#Right.
     *
     * @param k     Polynomial order (k ≥ 1). Sets K = k + 1 for dimensions.
     * @param t     Filter family/type code (e.g., Interpol, Legendre).
     *
     * @throws abort/error (via MRCPP messaging) on invalid @p k/@p t or if the
     *         required binary files cannot be opened.
     */
    CrossCorrelation(int k, int t);

    /**
     * @brief Construct from in-memory matrices (no file I/O).
     *
     * The order is inferred from the column count: 2K columns ⇒ order = K−1.
     * The two matrices must be shape-compatible (same size).
     *
     * @param t      Filter family/type code.
     * @param ldata  Left  matrix, size (K*K) × (2K).
     * @param rdata  Right matrix, size (K*K) × (2K).
     *
     * @throws abort/error if dimensions are inconsistent or the type is invalid.
     */
    CrossCorrelation(int t, const Eigen::MatrixXd &ldata, const Eigen::MatrixXd &rdata);

    /** @return The filter family/type code associated with this object. */
    int getType() const { return this->type; }

    /** @return The polynomial order k (so K = k + 1). */
    int getOrder() const { return this->order; }

    /** @return Const reference to the left cross-correlation matrix. */
    const Eigen::MatrixXd &getLMatrix() const { return this->Left; }

    /** @return Const reference to the right cross-correlation matrix. */
    const Eigen::MatrixXd &getRMatrix() const { return this->Right; }

protected:
    /**
     * @brief Filter family/type code.
     *
     * The meaning of this integer is validated against known families
     * (e.g., Interpolatory / Legendre) in the implementation. Kept as @c int
     * here to avoid header dependencies on the specific enum.
     */
    int type;

    /**
     * @brief Polynomial order k (k ≥ 1; K = k + 1).
     *
     * Controls the matrix dimensions:
     *   rows = K*K, cols = 2K.
     */
    int order;

    /**
     * @brief Left cross-correlation coefficient matrix.
     * Size: (K*K) × (2K), where K = order + 1.
     */
    Eigen::MatrixXd Left;

    /**
     * @brief Right cross-correlation coefficient matrix.
     * Size: (K*K) × (2K), where K = order + 1.
     */
    Eigen::MatrixXd Right;

private:
    /**
     * @brief Compose on-disk file paths for the left/right tables.
     *
     * Uses the discovered filter library root @p lib and the current
     * @c type / @c order to set #L_path and #R_path to the expected filenames.
     * (Naming convention is family-specific; see implementation.)
     */
    void setCCCPaths(const std::string &lib);

    /**
     * @brief Read the binary coefficient tables into #Left/#Right.
     *
     * Expects two files (left/right). Populates matrices with dimensions
     * (K*K) × (2K). Very small magnitudes may be zeroed for numerical
     * cleanliness (implementation detail).
     */
    void readCCCBin();

    /** @brief Full path to the left coefficient file (resolved at runtime). */
    std::string L_path;

    /** @brief Full path to the right coefficient file (resolved at runtime). */
    std::string R_path;
};

} // namespace mrcpp