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

#include <Eigen/Core>

#include "FunctionTree.h"
#include "MWNode.h"

namespace mrcpp {

template <int D, typename T> class FunctionNode final : public MWNode<D, T> {
public:
    FunctionTree<D, T> &getFuncTree() { return static_cast<FunctionTree<D, T> &>(*this->tree); }
    FunctionNode<D, T> &getFuncParent() { return static_cast<FunctionNode<D, T> &>(*this->parent); }
    FunctionNode<D, T> &getFuncChild(int i) { return static_cast<FunctionNode<D, T> &>(*this->children[i]); }

    const FunctionTree<D, T> &getFuncTree() const { return static_cast<const FunctionTree<D, T> &>(*this->tree); }
    const FunctionNode<D, T> &getFuncParent() const { return static_cast<const FunctionNode<D, T> &>(*this->parent); }
    const FunctionNode<D, T> &getFuncChild(int i) const { return static_cast<const FunctionNode<D, T> &>(*this->children[i]); }

    void createChildren(bool coefs) override;
    void genChildren() override;
    void genParent() override;
    void deleteChildren() override;

    T integrate() const;

    void setValues(const Eigen::Matrix<T , Eigen::Dynamic, 1> &vec);
    void getValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &vec);
    void getAbsCoefs(T *absCoefs);

    friend class FunctionTree<D, T>;
    friend class NodeAllocator<D, T>;

protected:
    FunctionNode()
            : MWNode<D, T>() {}
    FunctionNode(MWTree<D, T> *tree, int rIdx)
            : MWNode<D, T>(tree, rIdx) {}
    FunctionNode(MWNode<D, T> *parent, int cIdx)
            : MWNode<D, T>(parent, cIdx) {}
    FunctionNode(MWTree<D, T> *tree, const NodeIndex<D> &idx)
            : MWNode<D, T>(tree, idx) {}
    FunctionNode(const FunctionNode<D, T> &node) = delete;
    FunctionNode<D, T> &operator=(const FunctionNode<D, T> &node) = delete;
    ~FunctionNode() = default;

    T evalf(Coord<D> r);
    T evalScaling(const Coord<D> &r) const;

    void dealloc() override;
    void reCompress() override;

    T integrateLegendre() const;
    T integrateInterpolating() const;
    T integrateValues() const;
};

template <int D, typename T> T dot_scaling(const FunctionNode<D, T> &bra, const FunctionNode<D, T> &ket);
template <int D, typename T> T dot_wavelet(const FunctionNode<D, T> &bra, const FunctionNode<D, T> &ket);

} // namespace mrcpp
