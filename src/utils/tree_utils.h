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

#include "utils/math_utils.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {
namespace tree_utils {

template <int D, typename T> bool split_check(const MWNode<D, T> &node, double prec, double split_fac, bool abs_prec);

template <int D, typename T> void make_node_table(MWTree<D, T> &tree, MWNodeVector<D, T> &table);
template <int D, typename T> void make_node_table(MWTree<D, T> &tree, std::vector<MWNodeVector<D, T>> &table);

template <int D, typename T> void mw_transform(const MWTree<D, T> &tree, T *coeff_in, T *coeff_out, bool readOnlyScaling, int stride, bool overwrite = true);
//template <int D, typename T> void mw_transform_back(MWTree<D, T> &tree, T *coeff_in, T *coeff_out, int stride);
template <typename T> void mw_transform_back(MWTree<3, T> &tree, T *coeff_in, T *coeff_out, int stride);


} // namespace tree_utils
} // namespace mrcpp
