#pragma once

#include "mrcpp_declarations.h"

namespace mrcpp {
template<int D> void multiply(double prec, FunctionTree<D> &out, double c, FunctionTree<D> &tree_a, FunctionTree<D> &tree_b, int maxIter = -1);
template<int D> void multiply(double prec, FunctionTree<D> &out, FunctionTreeVector<D> &inp, int maxIter = -1);
template<int D> void power(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, double pow);
template<int D> void map(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, RepresentableFunction<D> &func);
}
