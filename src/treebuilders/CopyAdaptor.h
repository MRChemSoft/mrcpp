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

#include "TreeAdaptor.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {

template <int D> class CopyAdaptor final : public TreeAdaptor<D> {
public:
    CopyAdaptor(const FunctionTree<D> &t, int ms, int *bw);
    CopyAdaptor(const FunctionTreeVector<D> &t, int ms, int *bw);

private:
    int bandWidth[D];
    FunctionTreeVector<D> tree_vec;

    void setBandWidth(int *bw);
    bool splitNode(const MWNode<D> &node) const override;
};

} // namespace mrcpp
