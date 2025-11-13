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


/**
 * @brief Alias for a weighted FunctionTree pointer.
 *
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 *
 * @details
 * The tuple layout is:
 * - element 0: numeric coefficient of type @p T,
 * - element 1: pointer to a @c FunctionTree<D,T>.
 *
 * Ownership of the pointer is not implied by the alias; see @ref clear().
 */
template <int D, typename T = double>
using CoefsFunctionTree = std::tuple<T, FunctionTree<D, T> *>;

/**
 * @brief Alias for a vector of weighted FunctionTree pointers.
 *
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 */
template <int D, typename T = double>
using FunctionTreeVector = std::vector<CoefsFunctionTree<D, T>>;

/** @brief Remove all entries in the vector
 *  @param[in] fs: Vector to clear
 *  @param[in] dealloc: Option to free FunctionTree pointer before clearing
 */
template <int D, typename T>
void clear(FunctionTreeVector<D, T> &fs, bool dealloc = false) {
    if (dealloc) {
        for (auto &t : fs) {
            auto f = std::get<1>(t);
            if (f != nullptr) delete f;
            f = nullptr;
        }
    }
    fs.clear();
}

/**
 * @brief Compute the total number of nodes across all trees in the vector.
 *
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 * @param[in] fs Vector to fetch from
 * 
 * @returns Total number of nodes of all trees in the vector
 */
template <int D, typename T>
int get_n_nodes(const FunctionTreeVector<D, T> &fs) {
    int nNodes = 0;
    for (const auto &t : fs) {
        auto f = std::get<1>(t);
        if (f != nullptr) nNodes += f->getNNodes();
    }
    return nNodes;
}

/**
 * @brief Compute the total size of all trees in the vector (in kilobytes).
 * 
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 * @param[in] fs Vector to fetch from.
 * 
 * @returns Total size of all trees in the vector, in kB
 */
template <int D, typename T>
int get_size_nodes(const FunctionTreeVector<D, T> &fs) {
    int sNodes = 0;
    for (const auto &t : fs) {
        auto f = std::get<1>(t);
        if (f != nullptr) sNodes += f->getSizeNodes();
    }
    return sNodes;
}

/**
 * @brief Access the numeric coefficient at a given position.
 *
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 * @param[in] fs Vector to fetch from
 * @param[in] i  Position in vector (zero-based)
 * 
 * @returns Numerical coefficient at given position in vector
 *
 * @pre @p i must be a valid index in @p fs.
 */
template <int D, typename T>
T get_coef(const FunctionTreeVector<D, T> &fs, int i) {
    return std::get<0>(fs[i]);
}

/**
 * @brief Access the FunctionTree at a given position (non-const).
 *
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 * @param[in] fs Vector to fetch from
 * @param[in] i  Position in vector (zero-based)
 *
 * @return FunctionTree at given position in vector
 * 
 * @pre The pointer stored at position @p i must be non-null.
 */
template <int D, typename T>
FunctionTree<D, T> &get_func(FunctionTreeVector<D, T> &fs, int i) {
    return *(std::get<1>(fs[i]));
}

/**
 * @brief Access the FunctionTree at a given position (non-const).
 *
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 * @param[in] fs Vector to fetch from
 * @param[in] i  Position in vector (zero-based)
 *
 * @return FunctionTree at given position in vector
 * 
 * @pre The pointer stored at position @p i must be non-null.
 */
template <int D, typename T>
const FunctionTree<D, T> &get_func(const FunctionTreeVector<D, T> &fs, int i) {
    return *(std::get<1>(fs[i]));
}

} // namespace mrcpp