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

/**
 * @file CornerOperatorTree.h
 * @brief Declaration of CornerOperatorTree, a specialization of OperatorTree
 *        for "corner" non-standard form operators.
 *
 * @details
 * Many MRCPP operators are represented in non-standard form and decompose
 * into the four corner submatrices T, A, B, C. This helper class provides:
 * - computation of per-depth band widths for those corner blocks, and
 * - a fast band screen used during operator application.
 *
 * The band width information is stored in a BandWidth object owned by
 * the base class OperatorTree.
 *
 * @par Example
 * @code
 * CornerOperatorTree cot(mra, 10);            // maxDepth = 10
 * cot.calcBandWidth(1e-8);                    // build band widths with a threshold
 * bool within = cot.isOutsideBand(3, 4, 1);   // oTransl=3, o_depth=4, idx=1
 * @endcode
 */

/**
 * @class CornerOperatorTree
 * @brief Operator tree for non-standard form corner matrices.
 *
 * @details
 * This final class only adds band-handling logic on top of @ref OperatorTree.
 * Construction and storage are inherited from the base class; the only
 * public operations exposed here are:
 * - @ref calcBandWidth to build/update the band widths, and
 * - @ref isOutsideBand for a quick test against the stored band.
 */
class CornerOperatorTree final : public OperatorTree {
public:
    /// Inherit the valid constructor(s) from OperatorTree.
    using OperatorTree::OperatorTree;

    CornerOperatorTree(const CornerOperatorTree &tree) = delete;
    CornerOperatorTree &operator=(const CornerOperatorTree &tree) = delete;
    ~CornerOperatorTree() override = default;

    /**
     * @brief Compute per-depth band widths for the corner matrices.
     *
     * @param prec Threshold used when scanning matrix entries.
     *             If negative, the implementation falls back to the
     *             tree’s internal default (e.g. @c normPrec ).
     *
     * @details
     * For each depth and for each corner component \f$\{T,A,B,C\}\f$,
     * the routine scans along increasing distance from the diagonal and
     * records the largest translation \f$\ell\f$ for which the component
     * norm still exceeds the threshold. The resulting widths are stored
     * in the underlying @ref BandWidth structure.
     *
     * @note Calling this will (re)allocate and overwrite the stored band
     *       widths for the whole tree.
     */
    void calcBandWidth(double prec = -1.0) override;

    /**
     * @brief Test an offset against the stored band width.
     *
     * @param oTransl Integer offset (translation) from the diagonal.
     * @param o_depth Operator depth at which to query the band.
     * @param idx     Corner component selector in \f$\{0,1,2,3\}\f$
     *                corresponding to \f$\{T,A,B,C\}\f$.
     *
     * @return
     * **true** if \f$|oTransl| < \mathrm{width}(o\_depth, idx)\f$,
     * **false** otherwise.
     *
     * @details
     * Despite the historical name, this method returns @b true when the
     * offset lies @em inside the retained band (i.e., strictly smaller
     * than the stored width). Callers typically use it as a quick screen
     * to decide whether a sparse block needs to be applied.
     *
     * @warning This assumes @ref calcBandWidth has been called at least
     *          once; otherwise widths may be unset or conservative.
     */
    bool isOutsideBand(int oTransl, int o_depth, int idx) override;
};

} // namespace mrcpp