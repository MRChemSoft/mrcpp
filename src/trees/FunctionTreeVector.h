#pragma once

#include <tuple>
#include <vector>

#include "FunctionTree.h"

namespace mrcpp {

template <int D>
using CoefsFunctionTree = std::tuple<double, FunctionTree<D> *>;

template <int D>
using FunctionTreeVector = std::vector<CoefsFunctionTree<D>>;

template <int D>
void clear(FunctionTreeVector<D> & fs, bool dealloc = false) {
   if (dealloc) {
     for (auto & t : fs) {
       auto f = std::get<1>(t);
       if (f != nullptr) delete f;
     }
   }
   fs.clear();
}

template <int D>
int sumNodes(const FunctionTreeVector<D> & fs) {
  int nNodes = 0;
  for (const auto & t : fs) {
    auto f = std::get<1>(t);
    if (f != nullptr) {
      nNodes += f->getNNodes();
    }
  }
  return nNodes;
}

template <int D>
double getCoef(const FunctionTreeVector<D> & fs, int i) {
  return std::get<0>(fs[i]);
}

template <int D>
FunctionTree<D> & getFunc(FunctionTreeVector<D> & fs, int i) {
  return *(std::get<1>(fs[i]));
}

template <int D>
const FunctionTree<D> & getFunc(const FunctionTreeVector<D> & fs, int i) {
   return *(std::get<1>(fs[i]));
}
}
