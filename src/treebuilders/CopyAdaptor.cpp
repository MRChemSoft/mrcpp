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

#include "CopyAdaptor.h"

#include <tuple>

namespace mrcpp {

template <int D, typename T>
CopyAdaptor<D, T>::CopyAdaptor(FunctionTree<D, T> &t, int ms, int *bw)
        : TreeAdaptor<D, T>(ms) {
    setBandWidth(bw);
    tree_vec.push_back(std::make_tuple(1.0, &t));
}

template <int D, typename T>
CopyAdaptor<D, T>::CopyAdaptor(FunctionTreeVector<D, T> &t, int ms, int *bw)
        : TreeAdaptor<D, T>(ms)
        , tree_vec(t) {
    setBandWidth(bw);
}

template <int D, typename T> void CopyAdaptor<D, T>::setBandWidth(int *bw) {
    for (int d = 0; d < D; d++) {
        if (bw != nullptr) {
            this->bandWidth[d] = bw[d];
        } else {
            this->bandWidth[d] = 0;
        }
    }
}

template <int D, typename T> bool CopyAdaptor<D, T>::splitNode(const MWNode<D, T> &node) const {
    const NodeIndex<D> &idx = node.getNodeIndex();
    for (int c = 0; c < node.getTDim(); c++) {
        for (int d = 0; d < D; d++) {
            for (int bw = -this->bandWidth[d]; bw <= this->bandWidth[d]; bw++) {
                NodeIndex<D> bwIdx = idx.child(c);
                bwIdx[d] += bw;
                for (int i = 0; i < this->tree_vec.size(); i++) {
                    const FunctionTree<D, T> &func_i = get_func(tree_vec, i);
                    const MWNode<D, T> *node_i = func_i.findNode(bwIdx);
                    if (node_i != nullptr) return true;
                }
            }
        }
    }
    return false;
}

// Explicit instantiations
template class CopyAdaptor<1, double>;
template class CopyAdaptor<2, double>;
template class CopyAdaptor<3, double>;

template class CopyAdaptor<1, ComplexDouble>;
template class CopyAdaptor<2, ComplexDouble>;
template class CopyAdaptor<3, ComplexDouble>;

} // namespace mrcpp