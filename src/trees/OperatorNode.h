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
 * @file OperatorNode.h
 * @brief Node type for operator trees (2D non-standard form blocks).
 *
 * @details
 * This header declares @ref mrcpp::OperatorNode, a concrete node type used by
 * @ref mrcpp::OperatorTree. It specializes @ref mrcpp::MWNode with spatial
 * dimension @c D=2 to store and manage the coefficients of non-standard form
 * operator blocks (typically the corner blocks \f$T, A, B, C\f$).
 *
 * The class offers typed accessors to the owning @ref mrcpp::OperatorTree and
 * to parent/children nodes, and overrides a few hooks related to allocation,
 * child generation, and norm computation that are specific to operator nodes.
 */

#pragma once

#include <Eigen/Core>  // for Eigen::MatrixXd

#include "MWNode.h"
#include "OperatorTree.h"

namespace mrcpp {

/**
 * @class OperatorNode
 * @brief Leaf/branch node used inside @ref OperatorTree (fixed to 2D).
 *
 * @details
 * An OperatorNode stores the \f$2^D(k+1)^D\f$ multiwavelet coefficients for an
 * operator block at a given scale/translation and exposes:
 * - typed getters for the owning @ref OperatorTree and relatives,
 * - overrides for child creation/deletion and coefficient management,
 * - an overridden component-norm calculation suitable for operator blocks, and
 * - a helper to extract a single component block as an Eigen matrix.
 *
 * @note The spatial dimension is fixed to @c D=2 for operator trees.
 */
class OperatorNode final : public MWNode<2> {
public:
    /** @name Typed accessors to tree/relatives */
    ///@{
    /// @return Owning operator tree (non-const).
    OperatorTree &getOperTree() { return static_cast<OperatorTree &>(*this->tree); }
    /// @return Parent node as @ref OperatorNode (non-const).
    OperatorNode &getOperParent() { return static_cast<OperatorNode &>(*this->parent); }
    /// @return Child @p i as @ref OperatorNode (non-const).
    OperatorNode &getOperChild(int i) { return static_cast<OperatorNode &>(*this->children[i]); }

    /// @return Owning operator tree (const).
    const OperatorTree &getOperTree() const { return static_cast<const OperatorTree &>(*this->tree); }
    /// @return Parent node as @ref OperatorNode (const).
    const OperatorNode &getOperParent() const { return static_cast<const OperatorNode &>(*this->parent); }
    /// @return Child @p i as @ref OperatorNode (const).
    const OperatorNode &getOperChild(int i) const { return static_cast<const OperatorNode &>(*this->children[i]); }
    ///@}

    /**
     * @brief Create children nodes.
     * @param coefs If @c true, also allocate coefficient storage for each child.
     * @details Overrides @ref MWNode::createChildren to honor operator-specific
     * allocation and bookkeeping.
     */
    void createChildren(bool coefs) override;

    /**
     * @brief Generate children on demand (without necessarily allocating coefs).
     * @details Overrides @ref MWNode::genChildren to implement the operator-tree
     * generation policy.
     */
    void genChildren() override;

    /**
     * @brief Delete all children nodes (and their coefficient storage).
     * @details Overrides @ref MWNode::deleteChildren with operator-specific cleanup.
     */
    void deleteChildren() override;

    friend class OperatorTree;
    friend class NodeAllocator<2>;

protected:
    /** @name Construction and assignment */
    ///@{
    /// Default constructor (used by allocators).
    OperatorNode()
            : MWNode<2>(){};
    /// Root-node constructor (called by the owning tree).
    OperatorNode(MWTree<2> *tree, int rIdx)
            : MWNode<2>(tree, rIdx){};
    /// Child-node constructor (called when splitting a parent).
    OperatorNode(MWNode<2> *parent, int cIdx)
            : MWNode<2>(parent, cIdx){};
    /// Non-copyable.
    OperatorNode(const OperatorNode &node) = delete;
    /// Non-assignable.
    OperatorNode &operator=(const OperatorNode &node) = delete;
    /// Virtual destructor.
    ~OperatorNode() = default;
    ///@}

    /**
     * @brief Release coefficient storage (if owned) and reset node state.
     * @details Overrides @ref MWNode::dealloc to ensure operator-node invariants.
     */
    void dealloc() override;

    /**
     * @brief Compute squared norm of a specific component (scaling/wavelet block).
     * @param i Component index in \f$[0, 2^D)\f$.
     * @return Squared L2 norm of the requested component.
     * @details Overrides @ref MWNode::calcComponentNorm to match the operator
     * interpretation of components (e.g., corner blocks in non-standard form).
     */
    double calcComponentNorm(int i) const override;

    /**
     * @brief Extract a component block as a dense matrix.
     * @param i Component index in \f$[0, 2^D)\f$.
     * @return A matrix view/copy of the component coefficients (size \f$(k+1)\times(k+1)\f$).
     * @note Primarily intended for diagnostics and I/O.
     */
    Eigen::MatrixXd getComponent(int i);
};

} // namespace mrcpp