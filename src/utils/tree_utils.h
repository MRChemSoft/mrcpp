/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {
namespace tree_utils {

template<int D> bool split_check(const MWNode<D> &node, double prec, double split_fac, bool abs_prec);

template<int D> void make_node_table(MWTree<D> &tree, MWNodeVector<D> &table);
template<int D> void make_node_table(MWTree<D> &tree, std::vector<MWNodeVector<D>> &table);

template<int D> void mw_transform(const MWTree<D> &tree, double *coeff_in, double *coeff_out, bool readOnlyScaling, int stride, bool overwrite = true);
template<int D> void mw_transform_back(MWTree<D> &tree, double *coeff_in, double *coeff_out, int stride);

} // namespace node_utils
} // namespace mrcpp
