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

#include "MRCPP/mrcpp_declarations.h"
#include "trees/MWNode.h"

namespace mrcpp {

template <int D> class TreeAdaptor {
public:
    TreeAdaptor(int ms)
            : maxScale(ms) {}
    virtual ~TreeAdaptor() = default;

    void setMaxScale(int ms) { this->maxScale = ms; }

    void splitNodeVector(MWNodeVector<D> &out, MWNodeVector<D> &inp) const {
        for (int n = 0; n < inp.size(); n++) {
            MWNode<D> &node = *inp[n];
            // Can be BranchNode in operator application
            if (node.isBranchNode()) continue;
            if (node.getScale() + 2 > this->maxScale) continue;
            if (splitNode(node)) {
                node.createChildren(true);
                for (int i = 0; i < node.getNChildren(); i++) out.push_back(&node.getMWChild(i));
            }
        }
    }

protected:
    int maxScale;

    virtual bool splitNode(const MWNode<D> &node) const = 0;
};

} // namespace mrcpp
