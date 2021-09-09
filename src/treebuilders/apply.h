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

#include "trees/FunctionTreeVector.h"

namespace mrcpp {

// clang-format off
template <int D> class FunctionTree;
template <int D> class DerivativeOperator;
template <int D> class ConvolutionOperator;

template <int D> void apply(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter = -1, bool absPrec = false);
template <int D> void apply(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, FunctionTreeVector<D> &precTrees, int maxIter = -1, bool absPrec = false);
template <int D> void apply_far_field(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter = -1, bool absPrec = false);
template <int D> void apply_near_field(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter = -1, bool absPrec = false);
template <int D> void apply(FunctionTree<D> &out, DerivativeOperator<D> &oper, FunctionTree<D> &inp, int dir = -1);
template <int D> void divergence(FunctionTree<D> &out, DerivativeOperator<D> &oper, FunctionTreeVector<D> &inp);
template <int D> void divergence(FunctionTree<D> &out, DerivativeOperator<D> &oper, std::vector<FunctionTree<D> *> &inp);
template <int D> FunctionTreeVector<D> gradient(DerivativeOperator<D> &oper, FunctionTree<D> &inp);
// clang-format on

} // namespace mrcpp
