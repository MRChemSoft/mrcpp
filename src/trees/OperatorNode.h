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

template <typename T> class OperatorNodeT final : public MWNode<2, T> {
public:
    OperatorTreeT<T> &getOperTree() { return static_cast<OperatorTreeT<T> &>(*this->tree); }
    OperatorNodeT<T> &getOperParent() { return static_cast<OperatorNodeT<T> &>(*this->parent); }
    OperatorNodeT<T> &getOperChild(int i) { return static_cast<OperatorNodeT<T> &>(*this->children[i]); }

    const OperatorTreeT<T> &getOperTree() const { return static_cast<const OperatorTreeT<T> &>(*this->tree); }
    const OperatorNodeT<T> &getOperParent() const { return static_cast<const OperatorNodeT<T> &>(*this->parent); }
    const OperatorNodeT<T> &getOperChild(int i) const { return static_cast<const OperatorNodeT<T> &>(*this->children[i]); }

    void createChildren(bool coefs) override;
    void genChildren() override;
    void deleteChildren() override;

    friend class OperatorTreeT<T>;
    friend class NodeAllocator<2, T>;

protected:
    OperatorNodeT()
            : MWNode<2, T>(){};
    OperatorNodeT(MWTree<2, T> *tree, int rIdx)
            : MWNode<2, T>(tree, rIdx){};
    OperatorNodeT(MWNode<2, T> *parent, int cIdx)
            : MWNode<2, T>(parent, cIdx){};
    OperatorNodeT(const OperatorNodeT &node) = delete;
    OperatorNodeT &operator=(const OperatorNodeT &node) = delete;
    ~OperatorNodeT() = default;

    void dealloc() override;
    double calcComponentNorm(int i) const override;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getComponent(int i);
};

} // namespace mrcpp
