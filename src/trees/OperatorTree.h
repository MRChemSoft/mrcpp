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

#include <Eigen/Core>  // for Eigen::VectorXi

#include "MWTree.h"
#include "NodeAllocator.h"

namespace mrcpp {

// Forward declarations to avoid including the full headers here.
class BandWidth;
class OperatorNode;

/**
 * @class OperatorTree
 * @brief Base class for 2D operator trees in non-standard form
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
     * @brief Construct an operator tree
     * @param[in] mra  Multi-resolution analysis (domain + basis) shared by the tree
     * @param[in] np   “Norm precision” used when estimating/screening norms
     * @param[in] name Optional diagnostic name
     */
    OperatorTree(const MultiResolutionAnalysis<2> &mra, double np, const std::string &name = "nn");

    OperatorTree(const OperatorTree &tree) = delete;
    OperatorTree &operator=(const OperatorTree &tree) = delete;

    /// Virtual destructor
    virtual ~OperatorTree() override;

    /// @return The precision value used for norm-based screening 
    double getNormPrecision() const { return this->normPrec; }

    /**
     * @brief Release any existing @ref BandWidth object and set the pointer to null
     * @details Call this if the operator has changed and band widths must be recomputed
     */
    void clearBandWidth();

    /** @brief Calculates band widths of the non-standard form matrices
     *
     * @param[in] prec: Precision used for thresholding
     *
     * @details It is starting from \f$ l = 0 \f$ and updating the band width value each time we encounter
     * considerable value while keeping increasing \f$ l \f$, that stands for the distance to the diagonal
     */
    virtual void calcBandWidth(double prec = -1.0);

    /** @brief Checks if the distance to diagonal is bigger than the operator band width.
     *
     * @param[in] oTransl: distance to diagonal
     * @param[in] o_depth: scaling order
     * @param[in] idx: index corresponding to one of the matrices \f$ A, B, C \f$ or \f$ T \f$.ì
     *
     * @returns True if \b oTransl is outside of the band and False otherwise
     */
    virtual bool isOutsideBand(int oTransl, int o_depth, int idx);

    /** @brief Cleans up end nodes.
     *
     * @param[in] trust_scale: there is no cleaning down below \b trust_scale (it speeds up operator building).
     *
     * @details Traverses the tree and rewrites end nodes having branch node twins,
     * i. e. identical with respect to scale and translation.
     * This method is very handy, when an adaptive operator construction
     * can make a significunt noise at low scaling depth.
     * Its need comes from the fact that mwTransform up cannot override
     * rubbish that can potentially stick to end nodes at a particular level,
     * and as a result spread further up to the root with mwTransform.
     */
    void removeRoughScaleNoise(int trust_scale = 10);

    /**
     * @brief Make 1D lists, adressable from [-l, l] scale by scale, of operator node pointers for fast operator retrieval
     * 
     * @details Populates @ref nodePtrStore and @ref nodePtrAccess to avoid repeated lookups.
     * 
     * @warning This method is not thread safe,
     * since it projects missing operator nodes on the fly. Hence, it must NEVER
     * be called within a parallel region, or all hell will break loose. This is
     * not really a problem, but you have been warned.
     */
    void setupOperNodeCache();

    /// @brief Clear the operator-node caches built by @ref setupOperNodeCache()
    void clearOperNodeCache();

    /// @return Mutable reference to the stored @ref BandWidth (must exist)
    BandWidth &getBandWidth() { return *this->bandWidth; }
    /// @return Const reference to the stored @ref BandWidth (must exist)
    const BandWidth &getBandWidth() const { return *this->bandWidth; }

    /**
     * @brief Fast accessor to a node by indices (scale, diagonal distance)
     * 
     * @param[in] n Scale (depth measured from the root scale).
     * @param[in] l Distance to the diagonal (translation difference); l=0 hits the diagonal
     * 
     * @return Reference to the requested @ref OperatorNode.
     * @warning Valid only after calling @ref setupOperNodeCache()
     */
    OperatorNode &getNode(int n, int l) {
        return *nodePtrAccess[n][l];
    }
    /// @overload
    const OperatorNode &getNode(int n, int l) const { return *nodePtrAccess[n][l]; }






       
 
    /**
     * @brief Regenerate all s/d-coeffs by backtransformation, starting at the bottom and thus purifying all coefficients
     * 
     * @param overwrite If @c true, child coefficients may overwrite parent storage
     * 
     * @details Option to overwrite or add up existing
     * coefficients of BranchNodes (can be used after operator application).
     * Reimplementation of MWTree::mwTransform() without OMP, as calculation
     * of OperatorNorm is done using random vectors, which is non-deterministic
     * in parallel. FunctionTrees should be fine.
     */    
    void mwTransformUp() override;


    /**
     * @brief Regenerate all scaling coeffs by MW transformation of existing s/w-coeffs on coarser scales, starting at the rootNodes
     * 
     * @param overwrite If @c true, child coefficients may overwrite existing scaling coefficients
     * 
     * @details Option to overwrite or add up existing
     * coefficients of BranchNodes (can be used after operator application).
     * Reimplementation of MWTree::mwTransform() without OMP, as calculation
     * of OperatorNorm is done using random vectors, which is non-deterministic
     * in parallel. FunctionTrees should be fine.
     */
    void mwTransformDown(bool overwrite) override;



    
    /// @overload
    using MWTree<2>::getNode;
    /// @overload
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