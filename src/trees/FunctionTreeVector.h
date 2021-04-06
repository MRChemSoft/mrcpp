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

#include <tuple>
#include <vector>

#include "FunctionTree.h"

namespace mrcpp {

template <int D> using CoefsFunctionTree = std::tuple<double, FunctionTree<D> *>;
template <int D> using FunctionTreeVector = std::vector<CoefsFunctionTree<D>>;

/** @brief Remove all entries in the vector
 *  @param[in] fs: Vector to clear
 *  @param[in] dealloc: Option to free FunctionTree pointer before clearing
 */
template <int D> void clear(FunctionTreeVector<D> &fs, bool dealloc = false) {
    if (dealloc) {
        for (auto &t : fs) {
            auto f = std::get<1>(t);
            if (f != nullptr) delete f;
        }
    }
    fs.clear();
}

/** @returns Total number of nodes of all trees in the vector
 *  @param[in] fs: Vector to fetch from
 */
template <int D> int get_n_nodes(const FunctionTreeVector<D> &fs) {
    int nNodes = 0;
    for (const auto &t : fs) {
        auto f = std::get<1>(t);
        if (f != nullptr) nNodes += f->getNNodes();
    }
    return nNodes;
}

/** @returns Total size of all trees in the vector, in kB
 *  @param[in] fs: Vector to fetch from
 */
template <int D> int get_size_nodes(const FunctionTreeVector<D> &fs) {
    int sNodes = 0;
    for (const auto &t : fs) {
        auto f = std::get<1>(t);
        if (f != nullptr) sNodes += f->getSizeNodes();
    }
    return sNodes;
}

/** @returns Numerical coefficient at given position in vector
 *  @param[in] fs: Vector to fetch from
 *  @param[in] i: Position in vector
 */
template <int D> double get_coef(const FunctionTreeVector<D> &fs, int i) {
    return std::get<0>(fs[i]);
}

/** @returns FunctionTree at given position in vector
 *  @param[in] fs: Vector to fetch from
 *  @param[in] i: Position in vector
 */
template <int D> FunctionTree<D> &get_func(FunctionTreeVector<D> &fs, int i) {
    return *(std::get<1>(fs[i]));
}

/** @returns FunctionTree at given position in vector
 *  @param[in] fs: Vector to fetch from
 *  @param[in] i: Position in vector
 */
template <int D> const FunctionTree<D> &get_func(const FunctionTreeVector<D> &fs, int i) {
    return *(std::get<1>(fs[i]));
}
} // namespace mrcpp
