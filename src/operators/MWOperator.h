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

#include <vector>

#include "trees/MultiResolutionAnalysis.h"
#include "trees/OperatorTree.h"

namespace mrcpp {

/**
 * @class MWOperator
 * @brief Base class for multiwavelet (MW) operators with separated expansions.
 *
 * @tparam D Spatial dimension of the function space (1, 2, or 3).
 *
 * @details
 * An MW operator is represented as a (typically low-rank) separated expansion
 * whose per-term, per-dimension components are stored as pointers to
 * @ref OperatorTree objects. This class provides:
 *  - bookkeeping for the operator’s *root* scale and *reach* (bandwidth),
 *  - storage for the raw operator terms and their per-dimension assignments,
 *  - utilities for bandwidth analysis and component access, and
 *  - construction of the 2D operator-domain @ref MultiResolutionAnalysis used
 *    by @ref OperatorTree.
 *
 * Derived classes are responsible for building/populating @c raw_exp and then
 * calling @ref initOperExp() to map raw terms into directional components.
 */
template <int D> class MWOperator {
public:
    /**
     * @brief Construct an MW operator wrapper.
     *
     * @param mra   D-dimensional analysis describing the function space/domain.
     * @param root  Operator root level (coarsest level at which the operator lives).
     * @param reach Operator reach (half-width in levels at the root). Negative values
     *              can be interpreted by implementations as “auto”.
     */
    MWOperator(const MultiResolutionAnalysis<D> &mra, int root, int reach)
            : oper_root(root)
            , oper_reach(reach)
            , MRA(mra) {}

    MWOperator(const MWOperator &oper) = delete;            ///< Non-copyable
    MWOperator &operator=(const MWOperator &oper) = delete; ///< Non-assignable
    virtual ~MWOperator() = default;

    /**
     * @brief Number of separated terms currently active in the operator.
     */
    int size() const { return this->oper_exp.size(); }

    /**
     * @brief Maximum effective bandwidth at a given depth (scale).
     * @param depth Depth index; if negative, returns the maximum over all depths.
     * @return The maximum bandwidth, or -1 if @p depth is invalid.
     */
    int getMaxBandWidth(int depth = -1) const;

    /**
     * @brief Vector of maximum effective bandwidths per depth.
     * @return Reference to internal cache of maximum bandwidths.
     */
    const std::vector<int> &getMaxBandWidths() const { return this->band_max; }

    /**
     * @brief Compute effective bandwidths for all components at all depths.
     * @param prec Numeric tolerance used in bandwidth estimation.
     */
    void calcBandWidths(double prec);

    /**
     * @brief Clear cached bandwidth information in all components.
     */
    void clearBandWidths();

    /**
     * @brief Root level (coarsest scale) of the operator domain.
     */
    int getOperatorRoot() const { return this->oper_root; }

    /**
     * @brief Operator reach (half-width at the root level).
     */
    int getOperatorReach() const { return this->oper_reach; }

    /**
     * @brief Mutable access to the @p i-th separated term, @p d-th dimension component.
     * @param i Separated term index.
     * @param d Cartesian direction index (0..D-1).
     * @return Reference to the requested @ref OperatorTree.
     */
    OperatorTree &getComponent(int i, int d);

    /**
     * @brief Const access to the @p i-th separated term, @p d-th dimension component.
     * @param i Separated term index.
     * @param d Cartesian direction index (0..D-1).
     * @return Const reference to the requested @ref OperatorTree.
     */
    const OperatorTree &getComponent(int i, int d) const;

    /**
     * @brief Direct access to the array of D components for the @p i-th term.
     */
    std::array<OperatorTree *, D> &operator[](int i) { return this->oper_exp[i]; }

    /**
     * @brief Const direct access to the array of D components for the @p i-th term.
     */
    const std::array<OperatorTree *, D> &operator[](int i) const { return this->oper_exp[i]; }

protected:
    /** @name Operator geometry */
    ///@{
    int oper_root;   ///< Operator root level (coarsest scale).
    int oper_reach;  ///< Operator reach (half-width in levels at the root).
    MultiResolutionAnalysis<D> MRA; ///< Function-space analysis (domain and basis).
    ///@}

    /** @name Operator storage */
    ///@{
    std::vector<std::array<OperatorTree *, D>> oper_exp;      ///< Active separated terms by dimension.
    std::vector<std::unique_ptr<OperatorTree>> raw_exp;        ///< Owned raw operator terms (before assignment).
    std::vector<int> band_max;                                 ///< Maximum bandwidth per depth.
    ///@}

    /**
     * @brief Build the 2D operator-domain MRA used by @ref OperatorTree.
     * @details Operators act on a 2D lattice (row/column) even for D-D function spaces.
     */
    MultiResolutionAnalysis<2> getOperatorMRA() const;

    /**
     * @brief Initialize @ref oper_exp with @p M separated terms.
     * @details By default, assigns the first @p M raw terms isotropically across all D dimensions.
     * @param M Number of separated terms to activate.
     */
    void initOperExp(int M);

    /**
     * @brief Assign a particular operator component for term @p i and direction @p d.
     * @param i Term index in the separated expansion.
     * @param d Cartesian direction index (0..D-1).
     * @param oper Pointer to the @ref OperatorTree to be used for this component.
     */
    void assign(int i, int d, OperatorTree *oper) { this->oper_exp[i][d] = oper; }
};

} // namespace mrcpp