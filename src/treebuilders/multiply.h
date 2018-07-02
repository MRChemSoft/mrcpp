#pragma once

#include "trees/FunctionTreeVector.h"

namespace mrcpp {
template<int D> class RepresentableFunction;
template<int D> class FunctionTree;

template<int D> void multiply(double prec, FunctionTree<D> &out, double c, FunctionTree<D> &inp_a, FunctionTree<D> &inp_b, int maxIter = -1);
template<int D> void multiply(double prec, FunctionTree<D> &out, FunctionTreeVector<D> &inp, int maxIter = -1);
template<int D> void power(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, double pow, int maxIter = -1);
template<int D> void map(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, RepresentableFunction<D> &func);
template<int D> void dot(double prec, FunctionTree<D> &out, FunctionTreeVector<D> &inp_a, FunctionTreeVector<D> &inp_b, int maxIter = -1);
template<int D> double dot(FunctionTree<D> &bra, FunctionTree<D> &ket);

}
