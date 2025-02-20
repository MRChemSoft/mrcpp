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
#include "utils/CompFunction.h"

namespace mrcpp {

// clang-format off

template <int D, typename T> class FunctionTree;
template <int D> class DerivativeOperator;
template <int D> class ConvolutionOperator;

constexpr ComplexDouble  defaultMetric [4][4] ={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

template <int D, typename T> void apply(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, int maxIter = -1, bool absPrec = false);
template <int D> void apply(double prec, CompFunction<D> &out, ConvolutionOperator<D> &oper, const CompFunction<D> &inp, const ComplexDouble (*metric)[4] = defaultMetric, int maxIter = -1, bool absPrec = false);
template <int D, typename T> void apply(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, FunctionTreeVector<D, T> &precTrees, int maxIter = -1, bool absPrec = false);
template <int D, typename T> void apply(double prec, CompFunction<D> &out, ConvolutionOperator<D> &oper, CompFunction<D> &inp, FunctionTreeVector<D, T> *precTrees, ComplexDouble (*metric)[4] = nullptr, int maxIter = -1, bool absPrec = false);
template <int D, typename T> void apply_far_field(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, int maxIter = -1, bool absPrec = false);
template <int D> void apply_far_field(double prec, CompFunction<D> &out, ConvolutionOperator<D> &oper, CompFunction<D> &inp, const ComplexDouble (*metric)[4] = defaultMetric, int maxIter = -1, bool absPrec = false);
template <int D, typename T> void apply_near_field(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, int maxIter = -1, bool absPrec = false);
template <int D> void apply_near_field(double prec, CompFunction<D> &out, ConvolutionOperator<D> &oper, CompFunction<D> &inp, const ComplexDouble (*metric)[4] = defaultMetric, int maxIter = -1, bool absPrec = false);
template <int D, typename T> void apply(FunctionTree<D, T> &out, DerivativeOperator<D> &oper, FunctionTree<D, T> &inp, int dir = -1);
template <int D> void apply(CompFunction<D> &out, DerivativeOperator<D> &oper, CompFunction<D> &inp, int dir = -1, const ComplexDouble (*metric)[4] = defaultMetric);
template <int D, typename T> void divergence(FunctionTree<D, T> &out, DerivativeOperator<D> &oper, FunctionTreeVector<D, T> &inp);
template <int D, typename T> void divergence(CompFunction<D> &out, DerivativeOperator<D> &oper, FunctionTreeVector<D, T> *inp, const ComplexDouble (*metric)[4] = defaultMetric);
template <int D, typename T> void divergence(FunctionTree<D, T> &out, DerivativeOperator<D> &oper, std::vector<FunctionTree<D, T> *> &inp);
template <int D, typename T> void divergence(CompFunction<D> &out, DerivativeOperator<D> &oper, std::vector<FunctionTree<D, T> *> *inp, const ComplexDouble (*metric)[4]  = defaultMetric);
template <int D, typename T> FunctionTreeVector<D, T> gradient(DerivativeOperator<D> &oper, FunctionTree<D, T> &inp);
// template <int D>
std::vector<CompFunction<3>*> gradient(DerivativeOperator<3> &oper, CompFunction<3> &inp, ComplexDouble  (*metric)[4] = nullptr);
// clang-format on

} // namespace mrcpp
