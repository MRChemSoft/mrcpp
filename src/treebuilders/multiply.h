/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "trees/FunctionTreeVector.h"

namespace mrcpp {
template <int D> class RepresentableFunction;
template <int D> class FunctionTree;

template <int D> void dot(double prec,
                          FunctionTree<D> &out,
                          FunctionTreeVector<D> &inp_a,
                          FunctionTreeVector<D> &inp_b,
                          int maxIter = -1,
                          bool absPrec = false);

template <int D> double dot(FunctionTree<D> &bra,
                            FunctionTree<D> &ket);

template <int D> double node_norm_dot(FunctionTree<D> &bra,
                                      FunctionTree<D> &ket,
                                      bool exact = false);

template <int D> void multiply(double prec,
                               FunctionTree<D> &out,
                               double c,
                               FunctionTree<D> &inp_a,
                               FunctionTree<D> &inp_b,
                               int maxIter = -1,
                               bool absPrec = false,
                               bool useMaxNorms = false);

template <int D> void multiply(double prec,
                               FunctionTree<D> &out,
                               std::vector<FunctionTree<D> *> &inp,
                               int maxIter = -1,
                               bool absPrec = false,
                               bool useMaxNorms = false);

template <int D> void multiply(double prec,
                               FunctionTree<D> &out,
                               FunctionTreeVector<D> &inp,
                               int maxIter = -1,
                               bool absPrec = false,
                               bool useMaxNorms = false);

template <int D> void power(double prec,
                            FunctionTree<D> &out,
                            FunctionTree<D> &inp,
                            double p,
                            int maxIter = -1,
                            bool absPrec = false);

template <int D> void square(double prec,
                             FunctionTree<D> &out,
                             FunctionTree<D> &inp,
                             int maxIter = -1,
                             bool absPrec = false);

} // namespace mrcpp
