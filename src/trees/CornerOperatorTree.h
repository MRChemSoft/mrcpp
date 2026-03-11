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

#include "OperatorTree.h"

namespace mrcpp {

/** @class CornerOperatorTree
 *
 * @brief Special case of OperatorTree class
 *
 * @details Tree structure of operators having corner matrices
 * \f$ A, B, C \f$ in the non-standard form.
 */
class CornerOperatorTree final : public OperatorTree {
public:
    /// Inherit the valid constructorfrom OperatorTree.
    using OperatorTree::OperatorTree;

    CornerOperatorTree(const CornerOperatorTree &tree) = delete;
    CornerOperatorTree &operator=(const CornerOperatorTree &tree) = delete;
    ~CornerOperatorTree() override = default;

    /** 
     * @brief Calculates band widths of the non-standard form matrices
     * @param prec Precision used for thresholding
     *
     * @details It is starting from \f$ l = 2^n \f$ and updating the band width value each time we encounter
     * considerable value while keeping decreasing down to \f$ l = 0 \f$, that stands for the distance to the diagonal.
     * This procedure is repeated for each matrix \f$ A, B \f$ and \f$ C \f$.
     */
    void calcBandWidth(double prec = -1.0) override;

    /** 
     * @brief Checks if the distance to diagonal is lesser than the operator band width
     * @param oTransl distance to diagonal
     * @param o_depth scaling order
     * @param idx index corresponding to one of the matrices \f$ A, B, C \f$ or \f$ T \f$
     *
     * @returns True if @p oTransl is outside of the corner band (close to diagonal) and False otherwise.
     */
    bool isOutsideBand(int oTransl, int o_depth, int idx) override;
};

} // namespace mrcpp