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
 * @brief Container and loader for the left and right cross-correlation coefficient matrices
 *        of a multiwavelet filter family
 *
 * @details
 * Holds two dense matrices associated with the chosen polynomial order \f$ k \f$ (\f$ K = k+1 \f$) and
 * filter family #type (Interpol or Legendre):
 * \f[
 *   \mathrm{Left},\,\mathrm{Right} \in \mathbb{R}^{K^2 \times 2K}.
 * \f]
 * Each row corresponds to a flattened index pair \f$ (i,j) \f$ with \f$ i,j \in \{0,\ldots,K-1\} \f$
 * and each row stores a \f$ 2K \f$-wide correlation stencil. Objects are constructed either by loading
 * binary coefficient files from disk or by adopting in-memory matrices. There are no mutating public
 * methods; the class is a thread-safe value holder once constructed.
 *
 * @see CrossCorrelationCache for the process-wide singleton that caches these objects by order
 */
class CrossCorrelation final {
public:
    /**
     * @brief Construct by loading coefficient tables from the MRCPP filter library
     * @param k Polynomial order (\f$ k \geq 1 \f$; \f$ K = k+1 \f$)
     * @param t Filter family/type code (Interpol or Legendre)
     *
     * @details Discovers the library path via details::find_filters(), selects files by @p t and @p k,
     * and reads them into #Left and #Right
     *
     * @note Aborts via MSG_ABORT if @p k or @p t is invalid or the files cannot be opened
     */
    CrossCorrelation(int k, int t);

    /**
     * @brief Construct from in-memory matrices without file I/O
     * @param t     Filter family/type code
     * @param[in] ldata Left matrix of size \f$ K^2 \times 2K \f$
     * @param[in] rdata Right matrix of size \f$ K^2 \times 2K \f$ (must match dimensions of @p ldata)
     *
     * @details The polynomial order is inferred as \f$ k = \text{ldata.cols()}/2 - 1 \f$
     *
     * @note Aborts if dimensions are inconsistent or the type is invalid
     */
    CrossCorrelation(int t, const Eigen::MatrixXd &ldata, const Eigen::MatrixXd &rdata);

    /** @return The filter family/type code associated with this object */
    int getType() const { return this->type; }

    /** @return The polynomial order k (so K = k + 1) */
    int getOrder() const { return this->order; }

    /** @return Const reference to the left cross-correlation matrix */
    const Eigen::MatrixXd &getLMatrix() const { return this->Left; }

    /** @return Const reference to the right cross-correlation matrix */
    const Eigen::MatrixXd &getRMatrix() const { return this->Right; }

protected:
    int type;  ///< Filter family/type code (Interpol or Legendre; see constants.h)
    int order; ///< Polynomial order \f$ k \f$ (\f$ K = k+1 \f$; controls matrix dimensions \f$ K^2 \times 2K \f$)

    Eigen::MatrixXd Left;  ///< Left cross-correlation matrix of size \f$ K^2 \times 2K \f$
    Eigen::MatrixXd Right; ///< Right cross-correlation matrix of size \f$ K^2 \times 2K \f$

private:
    /**
     * @brief Compose the on-disk paths to the left and right binary tables
     * @param[in] lib Root directory of the MRCPP filter library
     *
     * @details Sets #L_path and #R_path using family-specific naming conventions based on #type and #order
     */
    void setCCCPaths(const std::string &lib);

    /**
     * @brief Read the binary coefficient tables from disk into #Left and #Right
     *
     * @details Reads two files whose paths are stored in #L_path and #R_path. Each file contains
     * \f$ K^2 \f$ rows of \f$ 2K \f$ doubles. Values smaller than MachinePrec are zeroed.
     */
    void readCCCBin();

    std::string L_path; ///< Full path to the left coefficient binary file (resolved at construction)
    std::string R_path; ///< Full path to the right coefficient binary file (resolved at construction)
};

} // namespace mrcpp