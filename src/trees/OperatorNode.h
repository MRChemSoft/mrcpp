/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

    void createChildren() override {
        MWNode<2>::createChildren();
        this->clearIsEndNode();
    }
    void genChildren() override {
        MWNode<2>::createChildren();
        this->clearIsEndNode();
        this->giveChildrenCoefs();
    }
    void deleteChildren() override {
        MWNode<2>::deleteChildren();
        this->setIsEndNode();
    }

    friend class MWNode<2>;
    friend class OperatorTree;
    friend class OperatorNodeAllocator;

protected:
    OperatorNode()
            : MWNode<2>(){};
    OperatorNode(const OperatorNode &node) = delete;
    OperatorNode &operator=(const OperatorNode &node) = delete;

    void dealloc() override;
    double calcComponentNorm(int i) const override;
};

} // namespace mrcpp
