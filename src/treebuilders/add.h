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

namespace mrcpp {

template <int D, typename T>
void add(double prec, FunctionTree<D, T> &out, T a, FunctionTree<D, T> &tree_a, T b, FunctionTree<D, T> &tree_b, int maxIter = -1, bool absPrec = false, bool conjugate = false);
template <int D, typename T> void add(double prec, FunctionTree<D, T> &out, FunctionTreeVector<D, T> &inp, int maxIter = -1, bool absPrec = false, bool conjugate = false);
template <int D, typename T> void add(double prec, FunctionTree<D, T> &out, std::vector<FunctionTree<D, T> *> &inp, int maxIter = -1, bool absPrec = false, bool conjugate = false);

} // namespace mrcpp
