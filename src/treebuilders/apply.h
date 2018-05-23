#pragma once

#include "trees/FunctionTreeVector.h"

namespace mrcpp {
template<int D> class FunctionTree;
template<int D> class DerivativeOperator;
template<int D> class ConvolutionOperator;

template<int D> void apply(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter = -1);
template<int D> void apply(FunctionTree<D> &out, DerivativeOperator<D> &oper, FunctionTree<D> &inp, int dir = -1);
template<int D> void divergence(FunctionTree<D> &out, DerivativeOperator<D> &oper, FunctionTreeVector<D> &inp);
template<int D> FunctionTreeVector<D> gradient(DerivativeOperator<D> &oper, FunctionTree<D> &inp);
}
