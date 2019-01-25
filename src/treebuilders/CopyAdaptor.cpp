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

#include "CopyAdaptor.h"

#include <tuple>

namespace mrcpp {

template <int D>
CopyAdaptor<D>::CopyAdaptor(FunctionTree<D> &t, int ms, int *bw)
        : TreeAdaptor<D>(ms) {
    setBandWidth(bw);
    tree_vec.push_back(std::make_tuple(1.0, &t));
}

template <int D>
CopyAdaptor<D>::CopyAdaptor(FunctionTreeVector<D> &t, int ms, int *bw)
        : TreeAdaptor<D>(ms)
        , tree_vec(t) {
    setBandWidth(bw);
}

template <int D> void CopyAdaptor<D>::setBandWidth(int *bw) {
    for (int d = 0; d < D; d++) {
        if (bw != nullptr) {
            this->bandWidth[d] = bw[d];
        } else {
            this->bandWidth[d] = 0;
        }
    }
}

template <int D> bool CopyAdaptor<D>::splitNode(const MWNode<D> &node) const {
    const NodeIndex<D> &idx = node.getNodeIndex();
    for (int c = 0; c < node.getTDim(); c++) {
        NodeIndex<D> cIdx(idx, c);
        for (int d = 0; d < D; d++) {
            for (int bw = -this->bandWidth[d]; bw <= this->bandWidth[d]; bw++) {
                NodeIndex<D> bwIdx(cIdx);
                bwIdx.getTranslation()[d] += bw;
                for (int i = 0; i < this->tree_vec.size(); i++) {
                    const FunctionTree<D> &func_i = get_func(tree_vec, i);
                    const MWNode<D> *node_i = func_i.findNode(bwIdx);
                    if (node_i != nullptr) return true;
                }
            }
        }
    }
    return false;
}

template class CopyAdaptor<1>;
template class CopyAdaptor<2>;
template class CopyAdaptor<3>;

} // namespace mrcpp
