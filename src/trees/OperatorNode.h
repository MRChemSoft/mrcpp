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

#include <Eigen/Core>  // for Eigen::MatrixXd

#include "MWNode.h"
#include "OperatorTree.h"

namespace mrcpp {

/**
 * @class OperatorNode
 *
 * @brief Node of an @ref OperatorTree.
 *
 * @details
 * An operator node is formally a 2D node which stores the coefficients of an operator
 * for a given scale and translation. The translation in this case corresponds to the difference 
 * in translation index betweeen input and output nodes of the function to which the operator is applied.
 * The scaling and wavelet structure of the nodes encodes the which component of the operator the coefficients
 * refer to (T scaling-scaling, C scaling-wavelet, B wavelet-scaling, A wavelet-wavelet) according to the non-standard form.
 *
 * The class offers typed accessors to the owning @ref mrcpp::OperatorTree and
 * to parent/children nodes, and overrides a few hooks related to allocation,
 * child generation, and norm computation that are specific to operator nodes.
 *
 * @note The spatial dimension is fixed to @c D=2 for operator trees.
 */
class OperatorNode final : public MWNode<2> {
public:
    OperatorTree &getOperTree() { return static_cast<OperatorTree &>(*this->tree); }                                ///< @return Owning operator tree    
    OperatorNode &getOperParent() { return static_cast<OperatorNode &>(*this->parent); }                            ///< @return Parent node
    OperatorNode &getOperChild(int i) { return static_cast<OperatorNode &>(*this->children[i]); }                   ///< @return Child @p i as @ref OperatorNode (non-const).
    const OperatorTree &getOperTree() const { return static_cast<const OperatorTree &>(*this->tree); }              ///< @return Owning operator tree
    const OperatorNode &getOperParent() const { return static_cast<const OperatorNode &>(*this->parent); }          ///< @return Parent node as @ref OperatorNode (const)
    const OperatorNode &getOperChild(int i) const { return static_cast<const OperatorNode &>(*this->children[i]); } ///< @return Child @p i as @ref OperatorNode (const).

    /**
     * @brief Create child nodes
     *
     * @param coefs If @c true, also allocate coefficient storage for each child.
     *
     * @details Overrides @ref MWNode::createChildren to honor operator-specific
     * allocation and bookkeeping.
     */
    void createChildren(bool coefs) override;

    /**
     * @brief Generate child nodes and populates their coefficients.
     *
     * @details Overrides @ref MWNode::genChildren to implement the operator
     */
    void genChildren() override;

    /**
     * @brief Delete all child nodes (and their coefficient storage)
     *
     * @details Overrides @ref MWNode::deleteChildren with operator-specific cleanup.
     */
    void deleteChildren() override;

    friend class OperatorTree;
    friend class NodeAllocator<2>;

protected:
    /** 
     * @brief Default constructor (used by allocators).
     */
    OperatorNode()
            : MWNode<2>(){};
    /** 
     * @brief Root node constructor (used by allocators).
     */
    OperatorNode(MWTree<2> *tree, int rIdx)
            : MWNode<2>(tree, rIdx){};
    /** 
     * @brief Child node constructor (used by allocators).
     */
    OperatorNode(MWNode<2> *parent, int cIdx)
            : MWNode<2>(parent, cIdx){};
    /** 
     * @brief Operator nodes cannot be copied.
     */
    OperatorNode(const OperatorNode &node) = delete;
    /** 
     * @brief Operator nodes cannot be assigned.
     */
    OperatorNode &operator=(const OperatorNode &node) = delete;
    /** 
     * @brief Default destructor.
     */
    ~OperatorNode() = default;
    /**
     * @brief Release coefficient storage (if owned) and reset node state.
     * @details Overrides @ref MWNode::dealloc to ensure operator-node invariants.
     */
    void dealloc() override;

    /**
     * @brief Calculate the norm of a given component of the OperatorNode
     *
     * @param[in] i: component index in [0, 3] (2D operator node has 4 components)
     *
     * @details OperatorNorms are defined as matrix 2-norms that are expensive to calculate.
     * Thus we calculate some cheaper upper bounds for this norm for thresholding.
     * First a simple vector norm, then a product of the 1- and infinity-norm.
     */
    double calcComponentNorm(int i) const override;

    /** @brief Gets a given component of the non-standard form.
     *
     * @param[in] i: Index enumerating the non-standard form component (A, B, C, T).
     * @returns The requested \f$ (k + 1) \times (k + 1) \f$-size matrix of the non-standard form.
     *
     * @details OperatorNode is uniquely associted with a scale \f$ n \f$ and translation
     * \f$ l = -2^n + 1, \ldots, 2^n = 1 \f$.
     * The non-standard form \f$ T_n, B_n, C_n, A_n \f$ defines matrices
     * \f$ \sigma_l^n, \beta_l^n, \gamma_l^n, \alpha_l^n \f$ for a given pair \f$ (n, l) \f$.
     * One of these matrices is returned by the method according to the choice of the index parameter
     * \f$ i = 0, 1, 2, 3 \f$, respectively.
     * For example, \f$ \alpha_l^n = \text{getComponent}(3) \f$.
     */
    Eigen::MatrixXd getComponent(int i);
};

} // namespace mrcpp
