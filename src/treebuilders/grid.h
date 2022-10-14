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

#include "functions/RepresentableFunction.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {
template <int D> void build_grid(FunctionTree<D> &out, int scales);
template <int D> void build_grid(FunctionTree<D> &out, const GaussExp<D> &inp, int maxIter = -1);
template <int D> void build_grid(FunctionTree<D> &out, const RepresentableFunction<D> &inp, int maxIter = -1);
template <int D> void build_grid(FunctionTree<D> &out, FunctionTree<D> &inp, int maxIter = -1);
template <int D> void build_grid(FunctionTree<D> &out, FunctionTreeVector<D> &inp, int maxIter = -1);
template <int D> void build_grid(FunctionTree<D> &out, std::vector<FunctionTree<D> *> &inp, int maxIter = -1);
template <int D> void copy_func(FunctionTree<D> &out, FunctionTree<D> &inp);
template <int D> void copy_grid(FunctionTree<D> &out, FunctionTree<D> &inp);
template <int D> void clear_grid(FunctionTree<D> &out);
template <int D> int refine_grid(FunctionTree<D> &out, int scales);
template <int D> int refine_grid(FunctionTree<D> &out, double prec, bool absPrec = false);
template <int D> int refine_grid(FunctionTree<D> &out, FunctionTree<D> &inp);
template <int D> int refine_grid(FunctionTree<D> &out, const RepresentableFunction<D> &inp);
} // namespace mrcpp
