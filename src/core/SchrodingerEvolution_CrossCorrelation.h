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
 * @class SchrodingerEvolution_CrossCorrelation
 * @brief Container and loader for cross-correlation coefficient matrices used
 *        by the Schrödinger time-evolution operator
 *
 * @details
 * Holds a sequence of \f$ (k+1) \times (k+1) \f$ matrices \f$ C^0, \ldots, C^{\text{amount}-1} \f$
 * encoding the scaling-function cross-correlation integrals needed to apply the
 * short-time Schrödinger propagator in the multiwavelet basis.
 * Coefficients are loaded from a binary file located via the MRCPP filter
 * library (details::find_filters()); currently only the Legendre scaling-basis
 * family is supported.
 *
 * @note This class is under active development. Member access specifiers
 *       (protected/private) are not yet finalized.
 *
 * @see CrossCorrelation for the conventional (non-time-evolution) counterpart
 */
class SchrodingerEvolution_CrossCorrelation final
{
public:
    /**
     * @brief Construct by loading coefficient matrices from the filter library
     *
     * @param amount The number of \f$ C^k \f$ matrices to load
     *               (\f$ k = 0, \ldots, \text{amount}-1 \f$)
     * @param k      Polynomial order of the scaling basis
     *               (must satisfy \f$ 1 \leq k \leq \f$ MaxOrder)
     * @param t      Scaling basis type tag (currently only Legendre is
     *               supported; Interpol is not yet implemented)
     *
     * @details
     * Validates @p k and @p t, resolves the binary coefficient file path
     * and reads the first @p amount matrices from disk. Each matrix is
     * truncated to the leading \f$ (k+1) \times (k+1) \f$ block.
     *
     * @throws Aborts via MSG_ABORT if @p k is out of range or the file
     *         cannot be opened. Emits MSG_ERROR for unsupported @p t values.
     */
    SchrodingerEvolution_CrossCorrelation(int amount, int k, int t);
    //SchrodingerEvolution_CrossCorrelation(int t, const Eigen::MatrixXd &ldata, const Eigen::MatrixXd &rdata);

    /** @return The scaling basis type code associated with this object */
    int getType() const { return this->type; }
    /** @return The polynomial order k (each matrix has size \f$ (k+1) \times (k+1) \f$) */
    int getOrder() const { return this->order; }
    /** @return The number of \f$ C^k \f$ matrices stored in #Matrix */
    int getAmount() const { return this->amount; }
    /** @return Const reference to the vector of stored \f$ C^k \f$ matrices */
    const std::vector<Eigen::MatrixXd> &getMatrix() const { return this->Matrix; }

//protected:
    int type;   ///< Scaling basis type code (e.g., Legendre)
    int order;  ///< Polynomial order k
    int amount; ///< Number of \f$ C^k \f$ matrices loaded and stored

    std::vector<Eigen::MatrixXd> Matrix; ///< Cross-correlation matrices \f$ C^0, \ldots, C^{\text{amount}-1} \f$

//private:
    /**
     * @brief Resolve the on-disk path for the coefficient binary
     * @param[in] lib Root directory of the MRCPP filter library
     *
     * Sets #path based on @c type; only the Legendre variant is currently
     * implemented
     */
    void setCCCPath(const std::string &lib);

    /**
     * @brief Read the binary file at #path and populate #Matrix
     *
     * Parses a leading Unicode text header (length-prefixed), followed by
     * the total matrix count and order, then reads the raw coefficient data.
     * Only the first #amount matrices are retained.
     */
    void readCCCBin();

    std::string path; ///< Full path to the binary cross-correlation coefficient file
};

} // namespace mrcpp
