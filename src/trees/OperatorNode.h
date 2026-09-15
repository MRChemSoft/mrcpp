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

#include "MWNode.h"
#include "OperatorTree.h"

namespace mrcpp {

/** @class OperatorNode
 *
 * @brief Node of an operator tree, holding one block of the non-standard form.
 *
 * @details Templated on the coefficient type alongside `OperatorTree`.
 */
template <typename T = double> class OperatorNode final : public MWNode<2, T> {
public:
    OperatorTree<T> &getOperTree() { return static_cast<OperatorTree<T> &>(*this->tree); }
    OperatorNode<T> &getOperParent() { return static_cast<OperatorNode<T> &>(*this->parent); }
    OperatorNode<T> &getOperChild(int i) { return static_cast<OperatorNode<T> &>(*this->children[i]); }

    const OperatorTree<T> &getOperTree() const { return static_cast<const OperatorTree<T> &>(*this->tree); }
    const OperatorNode<T> &getOperParent() const { return static_cast<const OperatorNode<T> &>(*this->parent); }
    const OperatorNode<T> &getOperChild(int i) const { return static_cast<const OperatorNode<T> &>(*this->children[i]); }

    void createChildren(bool coefs) override;
    void genChildren() override;
    void deleteChildren() override;

    /** @brief Matrix elements of the non-standard form.
     *
     * @param[in] i: index enumerating the matrix type in the non-standard form
     * @returns A \f$ (k + 1) \times (k + 1) \f$ submatrix of the non-standard form.
     *
     * @details An operator node is uniquely associated with a scale \f$ n \f$
     * and a translation \f$ l = -2^n + 1, \ldots, 2^n - 1 \f$. The
     * non-standard form defines matrices \f$ \sigma_l^n, \beta_l^n,
     * \gamma_l^n, \alpha_l^n \f$ for that pair, selected by
     * \f$ i = 0, 1, 2, 3 \f$ respectively.
     */
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getComponent(int i) const;

    friend class OperatorTree<T>;
    friend class NodeAllocator<2, T>;

protected:
    OperatorNode()
            : MWNode<2, T>(){};
    OperatorNode(MWTree<2, T> *tree, int rIdx)
            : MWNode<2, T>(tree, rIdx){};
    OperatorNode(MWNode<2, T> *parent, int cIdx)
            : MWNode<2, T>(parent, cIdx){};
    OperatorNode(const OperatorNode<T> &node) = delete;
    OperatorNode<T> &operator=(const OperatorNode<T> &node) = delete;
    ~OperatorNode() = default;

    void dealloc() override;
    double calcComponentNorm(int i) const override;
};

} // namespace mrcpp
