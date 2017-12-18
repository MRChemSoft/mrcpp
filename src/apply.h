#pragma once

#include "mrcpp_declarations.h"

namespace mrcpp {
template<int D> void apply(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter = -1);
template<int D> void apply(FunctionTree<D> &out, DerivativeOperator<D> &oper, FunctionTree<D> &inp, int dir = -1);
}
