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


#include "periodic_utils.h"

#include <cmath>

#include "details.h"
#include "Printer.h"
#include "trees/NodeIndex.h"

namespace mrcpp {
namespace periodic {

template <int D> bool in_unit_cell(const NodeIndex<D> &idx) {
    auto scale = idx.getScale();
    if (scale < 0) MSG_ABORT("Negative value in bit-shift");
    int two_n = 1 << scale;
    for (auto i = 0; i < D; i++) {
        if (idx[i] >= two_n) return false;
        if (idx[i] < 0) return false;
    }
    return true;
}

template <int D> NodeIndex<D> index_manipulation(const NodeIndex<D> &idx, const std::array<bool, D> &periodic) {
    if (!details::are_any(periodic, true)) return NodeIndex<D>(idx);

    auto scale = idx.getScale();
    if (scale < 0) {
        return NodeIndex<D>(scale);
    } else {
        int two_n = (1 << scale);
        std::array<int, D> new_l = idx.getTranslation();
        for (auto d = 0; d < D; ++d) {
            if (periodic[d]) {
                if (idx[d] >= two_n) new_l[d] = idx[d] % two_n;
                if (idx[d] < 0) new_l[d] = (idx[d] + 1) % two_n + two_n - 1;
            }
        }
        return NodeIndex<D>(scale, new_l);
    }
}

// Assumes a coordiate for a function defined on the [-1, 1] periodic cell.
// If r[i] is outside the unit-cell the function value is mapped to the unit-cell.
template <int D> Coord<D> coord_manipulation(const Coord<D> &r, const std::array<bool, D> &periodic) {
    if (!details::are_any(periodic, true)) return Coord<D>(r);

    Coord<D> new_r = r;
    for (auto d = 0; d < D; ++d) {
        if (periodic[d]) {
            if (r[d] >= 1.0) new_r[d] = std::fmod(r[d], 1.0);
            if (r[d] < 0.0) new_r[d] = std::fmod(r[d], 1.0) + 1.0;
        }
    }
    return new_r;
}

template bool in_unit_cell<1>(const NodeIndex<1> &idx);
template bool in_unit_cell<2>(const NodeIndex<2> &idx);
template bool in_unit_cell<3>(const NodeIndex<3> &idx);

template NodeIndex<1> index_manipulation<1>(const NodeIndex<1> &idx, const std::array<bool, 1> &periodic);
template NodeIndex<2> index_manipulation<2>(const NodeIndex<2> &idx, const std::array<bool, 2> &periodic);
template NodeIndex<3> index_manipulation<3>(const NodeIndex<3> &idx, const std::array<bool, 3> &periodic);

template Coord<1> coord_manipulation<1>(const Coord<1> &r, const std::array<bool, 1> &periodic);
template Coord<2> coord_manipulation<2>(const Coord<2> &r, const std::array<bool, 2> &periodic);
template Coord<3> coord_manipulation<3>(const Coord<3> &r, const std::array<bool, 3> &periodic);

} // namespace periodic
} // namespace mrcpp
