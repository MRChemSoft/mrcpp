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

template <int D> class FunctionNode final : public MWNode<D> {
public:
    FunctionTree<D> &getFuncTree() { return static_cast<FunctionTree<D> &>(*this->tree); }
    FunctionNode<D> &getFuncParent() { return static_cast<FunctionNode<D> &>(*this->parent); }
    FunctionNode<D> &getFuncChild(int i) { return static_cast<FunctionNode<D> &>(*this->children[i]); }

    const FunctionTree<D> &getFuncTree() const { return static_cast<const FunctionTree<D> &>(*this->tree); }
    const FunctionNode<D> &getFuncParent() const { return static_cast<const FunctionNode<D> &>(*this->parent); }
    const FunctionNode<D> &getFuncChild(int i) const { return static_cast<const FunctionNode<D> &>(*this->children[i]); }

    void createChildren(bool coefs) override;
    void genChildren() override;
    void genParent() override;
    void deleteChildren() override;

    double integrate() const;

    void setValues(const Eigen::VectorXd &vec);
    void getValues(Eigen::VectorXd &vec);
    void getAbsCoefs(double *absCoefs);

    friend class FunctionTree<D>;
    friend class NodeAllocator<D>;

protected:
    FunctionNode()
            : MWNode<D>() {}
    FunctionNode(MWTree<D> *tree, int rIdx)
            : MWNode<D>(tree, rIdx) {}
    FunctionNode(MWNode<D> *parent, int cIdx)
            : MWNode<D>(parent, cIdx) {}
    FunctionNode(MWTree<D> *tree, const NodeIndex<D> &idx)
            : MWNode<D>(tree, idx) {}
    FunctionNode(const FunctionNode<D> &node) = delete;
    FunctionNode<D> &operator=(const FunctionNode<D> &node) = delete;
    ~FunctionNode() = default;

    double evalf(Coord<D> r);
    double evalScaling(const Coord<D> &r) const;

    void dealloc() override;
    void reCompress() override;

    double integrateLegendre() const;
    double integrateInterpolating() const;
    double integrateValues() const;
};

template <int D> double dot_scaling(const FunctionNode<D> &bra, const FunctionNode<D> &ket);
template <int D> double dot_wavelet(const FunctionNode<D> &bra, const FunctionNode<D> &ket);

} // namespace mrcpp
