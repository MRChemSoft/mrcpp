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

#include "TreeCalculator.h"

namespace mrcpp {

template <int D, typename T> class SquareCalculator final : public TreeCalculator<D, T> {
public:
    SquareCalculator(FunctionTree<D, T> &inp)
            : func(&inp) {}

private:
    FunctionTree<D, T> *func;

    void calcNode(MWNode<D, T> &node_o) override {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        int n_coefs = node_o.getNCoefs();
        T *coefs_o = node_o.getCoefs();
        // This generates missing nodes
        MWNode<D, T> node_i = func->getNode(idx); // Copy node
        node_i.mwTransform(Reconstruction);
        node_i.cvTransform(Forward);
        const T *coefs_i = node_i.getCoefs();
        for (int j = 0; j < n_coefs; j++) { coefs_o[j] = coefs_i[j] * coefs_i[j]; }
        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

} // namespace mrcpp
