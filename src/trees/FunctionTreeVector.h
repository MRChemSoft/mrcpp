/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
int sum_nodes(const FunctionTreeVector<D> & fs) {
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
double get_coef(const FunctionTreeVector<D> & fs, int i) {
  return std::get<0>(fs[i]);
}

template <int D>
FunctionTree<D> & get_func(FunctionTreeVector<D> & fs, int i) {
  return *(std::get<1>(fs[i]));
}

template <int D>
const FunctionTree<D> & get_func(const FunctionTreeVector<D> & fs, int i) {
   return *(std::get<1>(fs[i]));
}
}
