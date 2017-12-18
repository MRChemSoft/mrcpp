#pragma once

#include "mrcpp_declarations.h"

namespace mrcpp {
template<int D> void add(double prec, FunctionTree<D> &out, double a, FunctionTree<D> &tree_a, double b, FunctionTree<D> &tree_b, int maxIter = -1);
template<int D> void add(double prec, FunctionTree<D> &out, FunctionTreeVector<D> &inp, int maxIter = -1);
}
