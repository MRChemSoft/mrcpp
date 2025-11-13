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
 * @class MWFilter
 * @brief Container for a 2K×2K multiwavelet filter bank and its block views.
 *
 * High-level
 * ----------
 * An MWFilter represents the matrix of a 1D multiwavelet transform for a given
 * polynomial order and family (type). With K = order + 1, the full transform
 * matrix has size 2K × 2K and is organized into four K × K blocks:
 *
 *     filter = [ G0  G1 ]   (top row: scaling channel)
 *              [ H0  H1 ]   (bottom row: wavelet channel)
 *
 * In the implementation (.cpp), G0/H0 are loaded from binary tables and
 * G1/H1 are derived by family-specific symmetry relations. Transposes of the
 * four blocks are also precomputed for the compression direction.
 *
 * Usage model
 * -----------
 * - Construct from (order, type) → loads data from disk and builds blocks.
 * - Construct from a given 2K×2K matrix → slices into blocks (no I/O).
 * - Multiply vectors/matrices with the transform or its transpose using
 *   apply()/applyInverse().
 * - Query individual K×K subfilters for compression or reconstruction.
 *
 * Notes on 'type'
 * ---------------
 * 'type' identifies the filter family (e.g., Interpol or Legendre). The exact
 * integer codes are defined elsewhere in MRCPP and validated in the .cpp.
 *
 * Dimension conventions
 * ---------------------
 * - order = k, K = k + 1
 * - Full transform: 2K × 2K (acts on 2K-length vectors / 2K-row matrices).
 */
class MWFilter final {
public:
    /**
     * @brief Construct from order and family type; loads blocks from disk.
     * @param k Polynomial order (k ≥ 0; with library-defined upper bound).
     * @param t Filter family/type tag (e.g., Interpol or Legendre).
     *
     * Side effects (see .cpp):
     *  - Locates binary tables on disk (family+order dependent).
     *  - Reads G0 and H0, synthesizes G1 and H1 by symmetry.
     *  - Assembles the full 2K×2K matrix 'filter'.
     */
    MWFilter(int k, int t);

    /**
     * @brief Construct directly from a full 2K×2K matrix (no I/O).
     * @param t    Filter family/type tag.
     * @param data Full transform matrix of size 2K×2K.
     *
     * The order is inferred as K = data.cols()/2, order = K - 1.
     * The four K×K blocks (and their transposes) are sliced from @p data.
     */
    MWFilter(int t, const Eigen::MatrixXd &data);

    /**
     * @name Apply the transform / its transpose
     * @{
     *
     * @brief Apply the forward/reconstruction transform: data ← filter * data.
     * Overloads exist for Eigen::MatrixXd and Eigen::VectorXd.
     *
     * @brief Apply the inverse/compression transform: data ← filter^T * data.
     * Overloads exist for Eigen::MatrixXd and Eigen::VectorXd.
     *
     * Precondition:
     *  - data.rows() must equal filter.cols() (i.e., 2K).
     */
    void apply(Eigen::MatrixXd &data) const;
    void apply(Eigen::VectorXd &data) const;
    void applyInverse(Eigen::MatrixXd &data) const;
    void applyInverse(Eigen::VectorXd &data) const;
    /** @} */

    /** @return Polynomial order k (so K = k + 1). */
    int getOrder() const { return this->order; }

    /** @return Filter family/type code. */
    int getType() const { return this->type; }

    /** @return Const reference to the full 2K×2K transform matrix. */
    const Eigen::MatrixXd &getFilter() const { return this->filter; }

    /**
     * @brief Return one of the four K×K subfilters.
     * @param i    Block index in the chosen operation's order (0..3).
     * @param oper Operation selector (direction), defaults to 0.
     *
     * Semantics (see .cpp):
     *  - For Reconstruction: blocks returned in order (H0, G0, H1, G1).
     *  - For Compression:    transposed blocks (H0^T, H1^T, G0^T, G1^T).
     *
     * The actual enum/integer values for 'oper' (e.g., Reconstruction/Compression)
     * are defined elsewhere (constants header). This method aborts on invalid
     * @p i or @p oper.
     */
    const Eigen::MatrixXd &getSubFilter(int i, int oper = 0) const;

    /**
     * @brief Shorthand: return the i-th compression subfilter (transposed form).
     *        Order: i=0→H0^T, 1→H1^T, 2→G0^T, 3→G1^T.
     */
    const Eigen::MatrixXd &getCompressionSubFilter(int i) const;

    /**
     * @brief Shorthand: return the i-th reconstruction subfilter (direct form).
     *        Order: i=0→H0, 1→G0, 2→H1, 3→G1.
     */
    const Eigen::MatrixXd &getReconstructionSubFilter(int i) const;

protected:
    /**
     * @brief Filter family/type tag (e.g., Interpol, Legendre).
     */
    int type;

    /**
     * @brief Polynomial order k (K = k + 1).
     */
    int order;

    /**
     * @brief Auxiliary dimension (reserved; may be unused in current code).
     */
    int dim;

    /**
     * @name Stored matrices
     * @{
     * @brief Full transform and its K×K sub-blocks (+ transposes).
     *
     * Layout: filter = [ G0  G1 ]
     *                 [ H0  H1 ]
     */
    Eigen::MatrixXd filter; ///< Full MW-transformation matrix
    Eigen::MatrixXd G0;
    Eigen::MatrixXd G1;
    Eigen::MatrixXd H0;
    Eigen::MatrixXd H1;
    // Transpose
    Eigen::MatrixXd G0t;
    Eigen::MatrixXd G1t;
    Eigen::MatrixXd H0t;
    Eigen::MatrixXd H1t;
    /** @} */

private:
    /**
     * @brief Compose on-disk paths to H0/G0 tables (family- and order-specific).
     *
     * Implemented in the .cpp; uses a discovered filter library root.
     * Sets #H_path and #G_path accordingly.
     */
    void setFilterPaths(const std::string &lib);

    /**
     * @brief Slice #filter into sub-blocks and compute their transposes.
     *
     * Used in the constructor that takes a full matrix.
     */
    void fillFilterBlocks();

    /**
     * @brief Load H0/G0 from disk, synthesize H1/G1 by symmetry, and
     *        precompute transposes. Populates #G0,#G1,#H0,#H1 and their ^T.
     */
    void generateBlocks();

    /** @brief Absolute file path to the H0 table. */
    std::string H_path;
    /** @brief Absolute file path to the G0 table. */
    std::string G_path;
};

} // namespace mrcpp