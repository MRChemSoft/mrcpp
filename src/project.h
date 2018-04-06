#pragma once

#include <functional>

#include "mrcpp_declarations.h"

namespace mrcpp {
template<int D> void project(double prec, FunctionTree<D> &out, RepresentableFunction<D> &inp, int maxIter = -1);
template<int D> void project(double prec, FunctionTree<D> &out, std::function<double (const double *r)> func, int maxIter = -1);
}
