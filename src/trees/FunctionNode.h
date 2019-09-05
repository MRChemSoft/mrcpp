/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

template <int D> class FunctionNode : public MWNode<D> {
public:
    FunctionTree<D> &getFuncTree() { return static_cast<FunctionTree<D> &>(*this->tree); }
    FunctionNode<D> &getFuncParent() { return static_cast<FunctionNode<D> &>(*this->parent); }
    FunctionNode<D> &getFuncChild(int i) { return static_cast<FunctionNode<D> &>(*this->children[i]); }

    const FunctionTree<D> &getFuncTree() const { return static_cast<const FunctionTree<D> &>(*this->tree); }
    const FunctionNode<D> &getFuncParent() const { return static_cast<const FunctionNode<D> &>(*this->parent); }
    const FunctionNode<D> &getFuncChild(int i) const {
        return static_cast<const FunctionNode<D> &>(*this->children[i]);
    }

    virtual void setValues(const Eigen::VectorXd &vec);
    virtual void getValues(Eigen::VectorXd &vec);
    virtual void getAbsCoefs(double *absCoefs);

    friend class FunctionTree<D>;

protected:
    FunctionNode()
            : MWNode<D>() {}
    FunctionNode(const FunctionNode<D> &node) = delete;
    FunctionTree<D> &operator=(const FunctionNode<D> &node) = delete;
    virtual ~FunctionNode() { assert(this->tree == 0); }

    double evalf(Coord<D> r);
    double evalScaling(const Coord<D> &r) const;

    double integrate() const;
    double integrateLegendre() const;
    double integrateInterpolating() const;
};

template <int D> double dotScaling(const FunctionNode<D> &bra, const FunctionNode<D> &ket);
template <int D> double dotWavelet(const FunctionNode<D> &bra, const FunctionNode<D> &ket);

} // namespace mrcpp
