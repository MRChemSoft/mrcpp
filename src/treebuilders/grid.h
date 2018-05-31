#pragma once

#include "functions/RepresentableFunction.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {
template<int D> void build_grid(FunctionTree<D> &out, const RepresentableFunction<D> &inp, int maxIter = -1);
template<int D> void build_grid(FunctionTree<D> &out, FunctionTree<D> &inp, int maxIter = -1);
template<int D> void build_grid(FunctionTree<D> &out, FunctionTreeVector<D> &inp, int maxIter = -1);
template<int D> void copy_func(FunctionTree<D> &out, FunctionTree<D> &inp);
template<int D> void copy_grid(FunctionTree<D> &out, FunctionTree<D> &inp);
template<int D> void clear_grid(FunctionTree<D> &out);
template<int D> int refine_grid(FunctionTree<D> &out, double prec);
template<int D> int refine_grid(FunctionTree<D> &out, FunctionTree<D> &inp);
}
