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

class OperatorNode final : public MWNode<2> {
public:
    OperatorTree &getOperTree() { return static_cast<OperatorTree &>(*this->tree); }
    OperatorNode &getOperParent() { return static_cast<OperatorNode &>(*this->parent); }
    OperatorNode &getOperChild(int i) { return static_cast<OperatorNode &>(*this->children[i]); }

    const OperatorTree &getOperTree() const { return static_cast<const OperatorTree &>(*this->tree); }
    const OperatorNode &getOperParent() const { return static_cast<const OperatorNode &>(*this->parent); }
    const OperatorNode &getOperChild(int i) const { return static_cast<const OperatorNode &>(*this->children[i]); }

    void createChildren(bool coefs) override;
    void genChildren() override;
    void deleteChildren() override;

    friend class OperatorTree;
    friend class NodeAllocator<2>;

protected:
    OperatorNode()
            : MWNode<2>(){};
    OperatorNode(MWTree<2> *tree, int rIdx)
            : MWNode<2>(tree, rIdx){};
    OperatorNode(MWNode<2> *parent, int cIdx)
            : MWNode<2>(parent, cIdx){};
    OperatorNode(const OperatorNode &node) = delete;
    OperatorNode &operator=(const OperatorNode &node) = delete;
    ~OperatorNode() = default;

    void dealloc() override;
public:
    double calcComponentNorm(int i) const override;
    Eigen::MatrixXd getComponent(int i);
};

} // namespace mrcpp
