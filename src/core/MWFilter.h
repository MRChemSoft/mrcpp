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
 * @brief Container for the \f$ 2K \times 2K \f$ multiwavelet filter matrix and its \f$ K \times K \f$ sub-blocks
 *
 * @details
 * An MWFilter represents the 1D multiwavelet transform matrix for a given polynomial order \f$ k \f$ and
 * scaling-basis family. With \f$ K = k+1 \f$, the full transform matrix is \f$ 2K \times 2K \f$ and is
 * partitioned into four \f$ K \times K \f$ blocks:
 * \f[
 *   \mathbf{F} = \begin{pmatrix} G_0 & G_1 \\ H_0 & H_1 \end{pmatrix}
 * \f]
 * where \f$ G_0, H_0 \f$ are read from binary tables on disk and \f$ G_1, H_1 \f$ are derived by
 * family-specific symmetry relations. Transposes of all four blocks are precomputed for the compression
 * (adjoint) direction. Construct from (order, type) to load data from disk, or from a pre-built
 * \f$ 2K \times 2K \f$ matrix to slice sub-blocks without file I/O.
 *
 * @see FilterCache for the process-wide singleton that caches MWFilter objects by order
 */
class MWFilter final {
public:
    /**
     * @brief Construct from polynomial order and filter family; loads sub-blocks from disk
     * @param k Polynomial order (\f$ k \geq 0 \f$, up to the library-defined MaxOrder)
     * @param t Filter family/type tag (Interpol or Legendre; validated at runtime)
     *
     * @details Locates binary tables on disk using the MRCPP filter library root, reads \f$ G_0 \f$ and
     * \f$ H_0 \f$, synthesizes \f$ G_1 \f$ and \f$ H_1 \f$ by family-specific symmetry, precomputes
     * transposed sub-blocks, and assembles the full \f$ 2K \times 2K \f$ matrix #filter.
     */
    MWFilter(int k, int t);

    /**
     * @brief Construct directly from a pre-built \f$ 2K \times 2K \f$ matrix without file I/O
     * @param t    Filter family/type tag
     * @param[in] data Full transform matrix of size \f$ 2K \times 2K \f$
     *
     * @details The polynomial order is inferred as \f$ k = \text{data.cols()}/2 - 1 \f$.
     * The four \f$ K \times K \f$ sub-blocks and their transposes are sliced from @p data.
     */
    MWFilter(int t, const Eigen::MatrixXd &data);

    /**
     * @brief Apply the reconstruction transform in-place: \f$ \mathbf{d} \leftarrow \mathbf{F}\,\mathbf{d} \f$
     * @param[in,out] data Matrix or vector of size \f$ 2K \times \cdot \f$ (rows must equal \f$ 2K \f$)
     */
    void apply(Eigen::MatrixXd &data) const;
    /** @brief Apply the reconstruction transform in-place (vector overload) */
    void apply(Eigen::VectorXd &data) const;
    /**
     * @brief Apply the compression (adjoint) transform in-place: \f$ \mathbf{d} \leftarrow \mathbf{F}^T\,\mathbf{d} \f$
     * @param[in,out] data Matrix or vector of size \f$ 2K \times \cdot \f$
     */
    void applyInverse(Eigen::MatrixXd &data) const;
    /** @brief Apply the compression (adjoint) transform in-place (vector overload) */
    void applyInverse(Eigen::VectorXd &data) const;

    /** @return Polynomial order k (so K = k + 1) */
    int getOrder() const { return this->order; }

    /** @return Filter family/type code */
    int getType() const { return this->type; }

    /** @return Const reference to the full 2K×2K transform matrix */
    const Eigen::MatrixXd &getFilter() const { return this->filter; }

    /**
     * @brief Return one of the four \f$ K \times K \f$ sub-blocks selected by operation and index
     * @param i    Sub-block index (0–3)
     * @param oper Operation selector: Reconstruction returns \f$(H_0,G_0,H_1,G_1)\f$ for
     *             \f$ i=0,1,2,3 \f$; Compression returns the transposed forms
     *             \f$(H_0^T,H_1^T,G_0^T,G_1^T)\f$
     * @return Const reference to the requested sub-block matrix
     *
     * @note Aborts if @p i or @p oper is out of range
     */
    const Eigen::MatrixXd &getSubFilter(int i, int oper = 0) const;

    /**
     * @brief Return the @p i-th compression sub-block (transposed form)
     * @param i Index in \f$\{0,1,2,3\}\f$ selecting \f$H_0^T, H_1^T, G_0^T, G_1^T\f$ respectively
     * @return Const reference to the selected transposed sub-block
     */
    const Eigen::MatrixXd &getCompressionSubFilter(int i) const;

    /**
     * @brief Return the @p i-th reconstruction sub-block (direct form)
     * @param i Index in \f$\{0,1,2,3\}\f$ selecting \f$H_0, G_0, H_1, G_1\f$ respectively
     * @return Const reference to the selected sub-block
     */
    const Eigen::MatrixXd &getReconstructionSubFilter(int i) const;

protected:
    int type;  ///< Filter family/type tag (Interpol or Legendre; see constants.h)
    int order; ///< Polynomial order \f$ k \f$ (so \f$ K = k+1 \f$)
    int dim;   ///< Auxiliary dimension (reserved; currently unused)

    Eigen::MatrixXd filter; ///< Full \f$ 2K \times 2K \f$ MW transform matrix
    Eigen::MatrixXd G0;     ///< Scaling sub-block (top-left, size \f$ K \times K \f$)
    Eigen::MatrixXd G1;     ///< Scaling sub-block (top-right, size \f$ K \times K \f$)
    Eigen::MatrixXd H0;     ///< Wavelet sub-block (bottom-left, size \f$ K \times K \f$)
    Eigen::MatrixXd H1;     ///< Wavelet sub-block (bottom-right, size \f$ K \times K \f$)
    Eigen::MatrixXd G0t;    ///< Transpose of \f$ G_0 \f$
    Eigen::MatrixXd G1t;    ///< Transpose of \f$ G_1 \f$
    Eigen::MatrixXd H0t;    ///< Transpose of \f$ H_0 \f$
    Eigen::MatrixXd H1t;    ///< Transpose of \f$ H_1 \f$

private:
    /**
     * @brief Compose the on-disk paths to the \f$ H_0 \f$ and \f$ G_0 \f$ binary tables
     * @param[in] lib Root directory of the MRCPP filter library
     *
     * @details Sets #H_path and #G_path based on #type and #order using family-specific naming conventions
     */
    void setFilterPaths(const std::string &lib);

    /**
     * @brief Slice #filter into the four \f$ K \times K \f$ sub-blocks and compute their transposes
     *
     * @details Used by the matrix-constructor to populate #G0, #G1, #H0, #H1 and their transposed counterparts
     */
    void fillFilterBlocks();

    /**
     * @brief Load \f$ G_0 \f$ and \f$ H_0 \f$ from disk, synthesize \f$ G_1 \f$ and \f$ H_1 \f$ by symmetry,
     * and precompute all transposed sub-blocks
     *
     * @details Family-specific symmetry relations (Interpol vs. Legendre) determine how \f$ G_1 \f$ and
     * \f$ H_1 \f$ are derived from the loaded sub-blocks
     */
    void generateBlocks();

    std::string H_path; ///< Absolute path to the binary \f$ H_0 \f$ table file
    std::string G_path; ///< Absolute path to the binary \f$ G_0 \f$ table file
};

} // namespace mrcpp