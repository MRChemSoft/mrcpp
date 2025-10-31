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
 * @file OperatorTree.h
 * @brief Declaration of the multiwavelet operator tree (2D non-standard form).
 *
 * @details
 * An @ref mrcpp::OperatorTree stores a bivariate (D=2) operator in
 * multiwavelet (MW) **non-standard form**, i.e. split into corner blocks
 * \f$T, A, B, C\f$ at each scale. It provides:
 * - adaptive storage and traversal via the base @ref mrcpp::MWTree,
 * - optional **band-width screening** of corner blocks through @ref BandWidth,
 * - cached direct access to operator nodes to avoid repeated tree lookups, and
 * - MW up/down transforms specialized for operator data.
 *
 * Only trees built from **compatible** MRAs (same domain, order, and depth)
 * should be combined in further computations.
 */

#pragma once

#include <Eigen/Core>  // for Eigen::VectorXi

#include "MWTree.h"
#include "NodeAllocator.h"

namespace mrcpp {

// Forward declarations to avoid including the full headers here.
class BandWidth;
class OperatorNode;

/**
 * @class OperatorTree
 * @brief Base class for 2D operator trees in non-standard form.
 *
 * @details
 * The tree is organized like any MW tree (roots/branches/leaves) but stores
 * operator coefficients. A per-depth **band width** (distance from the main
 * diagonal in translation space) can be estimated to prune negligible corner
 * blocks during application.
 */
class OperatorTree : public MWTree<2> {
public:
    /**
     * @brief Construct an operator tree.
     * @param mra  Multi-resolution analysis (domain + basis) shared by the tree.
     * @param np   “Norm precision” used when estimating/screening norms.
     * @param name Optional diagnostic name.
     */
    OperatorTree(const MultiResolutionAnalysis<2> &mra, double np, const std::string &name = "nn");

    OperatorTree(const OperatorTree &tree) = delete;
    OperatorTree &operator=(const OperatorTree &tree) = delete;

    /// Virtual destructor.
    virtual ~OperatorTree() override;

    /// @return The precision value used for norm-based screening.
    double getNormPrecision() const { return this->normPrec; }

    /**
     * @brief Release any existing @ref BandWidth object and set the pointer to null.
     * @details Call this if the operator has changed and band widths must be recomputed.
     */
    void clearBandWidth();

    /**
     * @brief Estimate per-depth band widths for the corner matrices.
     * @param prec Threshold used when deciding if a component is significant.
     *             If negative, the internal @ref getNormPrecision() is used.
     * @details Populates the internally owned @ref BandWidth structure.
     */
    virtual void calcBandWidth(double prec = -1.0);

    /**
     * @brief Quick band-screening predicate.
     * @param oTransl  Distance from the diagonal in translation space (|l\_bra−l\_ket|).
     * @param o_depth  Depth/scale index where the test is performed.
     * @param idx      Corner block selector: 0 = T, 1 = C, 2 = B, 3 = A (convention as used internally).
     * @return @c true if @p oTransl is **outside** the currently stored band at @p o_depth for block @p idx.
     * @note Requires a previously computed @ref BandWidth (see @ref calcBandWidth()).
     */
    virtual bool isOutsideBand(int oTransl, int o_depth, int idx);

    /**
     * @brief Dampen/remove rough-scale numerical noise in the operator.
     * @param trust_scale Scales finer (greater or equal to this) are trusted and preserved.
     * @details Useful after building operators from noisy input data.
     */
    void removeRoughScaleNoise(int trust_scale = 10);

    /**
     * @brief Build cache tables for direct @ref OperatorNode access.
     * @details Populates @ref nodePtrStore and @ref nodePtrAccess to avoid repeated lookups.
     */
    void setupOperNodeCache();

    /// @brief Clear the operator-node caches built by @ref setupOperNodeCache().
    void clearOperNodeCache();

    /// @return Mutable reference to the stored @ref BandWidth (must exist).
    BandWidth &getBandWidth() { return *this->bandWidth; }
    /// @return Const reference to the stored @ref BandWidth (must exist).
    const BandWidth &getBandWidth() const { return *this->bandWidth; }

    /**
     * @brief Fast accessor to a node by (scale, diagonal distance).
     * @param n Scale (depth measured from the root scale).
     * @param l Distance to the diagonal (translation difference); l=0 hits the diagonal.
     * @return Reference to the requested @ref OperatorNode.
     * @warning Valid only after calling @ref setupOperNodeCache().
     */
    OperatorNode &getNode(int n, int l) {
        return *nodePtrAccess[n][l];
    }
    /// Const overload of @ref getNode(int,int).
    const OperatorNode &getNode(int n, int l) const { return *nodePtrAccess[n][l]; }

    /**
     * @brief Downward MW transform specialized for operator data.
     * @param overwrite If @c true, child coefficients may overwrite parent storage.
     */
    void mwTransformDown(bool overwrite) override;

    /// @brief Upward MW transform specialized for operator data.
    void mwTransformUp() override;

    // Bring MWTree overloads into scope.
    using MWTree<2>::getNode;
    using MWTree<2>::findNode;

protected:
    const double normPrec;     ///< Default precision used in norm-based heuristics.
    BandWidth *bandWidth;      ///< Optional per-depth band-width description (owned).

    /// @name Operator-node cache (built by @ref setupOperNodeCache()).
    ///@{
    OperatorNode ***nodePtrStore;   ///< Storage for contiguous (n,l) -> node pointers.
    OperatorNode ***nodePtrAccess;  ///< Centered view so that index l=0 addresses the diagonal.
    ///@}

    /// @brief Allocate all root nodes required by the current world box.
    void allocRootNodes();

    /**
     * @brief Compute the maximum translation index at each depth.
     * @param[out] maxTransl Vector whose @c d-th entry stores the maximum |l| at that depth.
     */
    void getMaxTranslations(Eigen::VectorXi &maxTransl);

    /// @brief Human-readable dump of tree statistics.
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp