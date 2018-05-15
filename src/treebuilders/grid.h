#pragma once

#include "functions/RepresentableFunction.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {
template<int D> void build_grid(FunctionTree<D> &out, const RepresentableFunction<D> &inp, int maxIter = -1);
template<int D> void copy_grid(FunctionTree<D> &out, FunctionTree<D> &inp, int maxIter = -1);
template<int D> void copy_grid(FunctionTree<D> &out, FunctionTreeVector<D> &inp, int maxIter = -1);
template<int D> int clear_grid(double prec, FunctionTree<D> &out);
}
