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

/**
 * @file MWOperator.cpp
 * @brief Common utilities for multiwavelet (MW) operators: term assignment,
 *        component access, bandwidth analysis, and operator-domain MRA setup.
 *
 * @details
 * The templated @ref mrcpp::MWOperator provides infrastructure shared by concrete
 * MW operators:
 *  - organizing a separated operator expansion into per-dimension @ref OperatorTree
 *    components,
 *  - computing effective bandwidths across scales, and
 *  - constructing the 2D operator-domain @ref MultiResolutionAnalysis used by
 *    operator trees (for D-dimensional function spaces).
 *
 * The functions implemented here are thin, performance-oriented utilities that
 * avoid modifying operator semantics. They are used by higher-level operators
 * such as convolution- and derivative-based classes.
 */

#include "MWOperator.h"
#include "trees/BandWidth.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

using namespace Eigen;

namespace mrcpp {

/**
 * @brief Initialize the separated operator expansion with @p M terms.
 *
 * @tparam D Spatial dimension of the target space.
 * @param M  Number of separated terms to activate from @c raw_exp.
 *
 * @details
 * Allocates an @c oper_exp array of size @p M, each entry holding @c D pointers
 * to @ref OperatorTree components (one per Cartesian direction).
 * By default, an *isotropic* operator is formed by assigning the first @p M raw
 * terms to **all** directions.
 *
 * @pre @c raw_exp has at least @p M terms; otherwise the function aborts.
 * @note You can override individual components later via @ref assign().
 */
template <int D> void MWOperator<D>::initOperExp(int M) {
    if (this->raw_exp.size() < M) MSG_ABORT("Incompatible raw expansion");
    this->oper_exp.clear();
    for (int m = 0; m < M; m++) {
        std::array<OperatorTree *, D> otrees;
        otrees.fill(nullptr);
        this->oper_exp.push_back(otrees);
    }

    // Sets up an isotropic operator with the first M raw terms in all directions
    for (int i = 0; i < M; i++)
        for (int d = 0; d < D; d++) assign(i, d, this->raw_exp[i].get());
}

/**
 * @brief Mutable access to a specific separated component.
 *
 * @tparam D Spatial dimension.
 * @param i  Term index in the separated expansion.
 * @param d  Cartesian direction index (0..D-1).
 * @return Reference to the requested @ref OperatorTree.
 *
 * @throws If indices are out of bounds or the component is null.
 */
template <int D> OperatorTree &MWOperator<D>::getComponent(int i, int d) {
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Index out of bounds");
    if (d < 0 or d >= D) MSG_ERROR("Dimension out of bounds");
    if (this->oper_exp[i][d] == nullptr) MSG_ERROR("Invalid component");
    return *this->oper_exp[i][d];
}

/**
 * @brief Const access to a specific separated component.
 *
 * @tparam D Spatial dimension.
 * @param i  Term index in the separated expansion.
 * @param d  Cartesian direction index (0..D-1).
 * @return Const reference to the requested @ref OperatorTree.
 *
 * @throws If indices are out of bounds or the component is null.
 */
template <int D> const OperatorTree &MWOperator<D>::getComponent(int i, int d) const {
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Index out of bounds");
    if (d < 0 or d >= D) MSG_ERROR("Dimension out of bounds");
    if (this->oper_exp[i][d] == nullptr) MSG_ERROR("Invalid component");
    return *this->oper_exp[i][d];
}

/**
 * @brief Get the maximum effective bandwidth at a given depth.
 *
 * @tparam D Spatial dimension.
 * @param depth Tree depth (scale) at which to query the bandwidth. If negative,
 *              the maximum over all depths is returned.
 * @return Non-negative maximum bandwidth, or -1 if @p depth is invalid.
 */
template <int D> int MWOperator<D>::getMaxBandWidth(int depth) const {
    int maxWidth = -1;
    if (depth < 0) {
        maxWidth = *std::max_element(this->band_max.begin(), this->band_max.end());
    } else if (depth < this->band_max.size()) {
        maxWidth = this->band_max[depth];
    }
    return maxWidth;
}

/**
 * @brief Clear cached @ref BandWidth information in all operator components.
 *
 * @tparam D Spatial dimension.
 */
template <int D> void MWOperator<D>::clearBandWidths() {
    for (auto &i : this->oper_exp)
        for (int d = 0; d < D; d++) i[d]->clearBandWidth();
}

/**
 * @brief Compute effective bandwidths at all scales for all components.
 *
 * @tparam D Spatial dimension.
 * @param prec Numerical precision used to estimate bandwidths.
 *
 * @details
 * For each @ref OperatorTree component, this calls @ref OperatorTree::calcBandWidth
 * and records the @em maximum effective width across components for every depth.
 * Results are stored in @c band_max and summarized to the log at print level 20.
 */
template <int D> void MWOperator<D>::calcBandWidths(double prec) {
    int maxDepth = 0;
    // First compute BandWidths and find depth of the deepest component
    for (auto &i : this->oper_exp) {
        for (int d = 0; d < D; d++) {
            OperatorTree &oTree = *i[d];
            oTree.calcBandWidth(prec);
            const BandWidth &bw = oTree.getBandWidth();
            int depth = bw.getDepth();
            if (depth > maxDepth) maxDepth = depth;
        }
    }
    this->band_max = std::vector<int>(maxDepth + 1, -1);

    // Find the largest effective bandwidth at each scale
    for (auto &i : this->oper_exp) {
        for (int d = 0; d < D; d++) {
            const OperatorTree &oTree = *i[d];
            const BandWidth &bw = oTree.getBandWidth();
            for (int n = 0; n <= bw.getDepth(); n++) { // scale loop
                for (int j = 0; j < 4; j++) {          // component loop
                    int w = bw.getWidth(n, j);
                    if (w > this->band_max[n]) this->band_max[n] = w;
                }
            }
        }
    }
    println(20, "  Maximum bandwidths:");
    for (auto bw : this->band_max) println(20, bw);
    println(20, std::endl);
}

/**
 * @brief Build the 2D operator-domain MRA used by operator trees.
 *
 * @tparam D Spatial dimension of the *function* domain.
 * @return A @ref MultiResolutionAnalysis<2> describing the operator lattice.
 *
 * @details
 * Operator trees live on a 2D lattice (row/column), even when acting on
 * D-dimensional function spaces. The lattice extents are determined from the
 * operator's root level and reach, and it reuses the function-space scaling
 * basis (uniform scaling is assumed).
 */
template <int D> MultiResolutionAnalysis<2> MWOperator<D>::getOperatorMRA() const {
    const BoundingBox<D> &box = this->MRA.getWorldBox();
    const ScalingBasis &basis = this->MRA.getScalingBasis();

    int reach = this->oper_reach + 1;
    if (reach < 0) {
        for (int i = 0; i < D; i++) {
            if (box.size(i) > reach) reach = box.size(i);
        }
    }
    auto l = std::array<int, 2>{};
    auto nbox = std::array<int, 2>{reach, reach};
    // Zero in argument since operators are only implemented
    // for uniform scaling factor
    auto sf = std::array<double, 2>{box.getScalingFactor(0), box.getScalingFactor(0)};

    BoundingBox<2> oper_box(this->oper_root, l, nbox, sf);
    auto oper_mra = MultiResolutionAnalysis<2>(oper_box, basis);
    return oper_mra;
}

/* Explicit template instantiations */
template class MWOperator<1>;
template class MWOperator<2>;
template class MWOperator<3>;

} // namespace mrcpp