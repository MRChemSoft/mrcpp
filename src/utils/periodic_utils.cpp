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

#include "MRCPP/Printer"

#include "periodic_utils.h"
#include "trees/NodeIndex.h"

#include <cmath>

namespace mrcpp {
namespace periodic {

template <int D> bool in_unit_cell(NodeIndex<D> idx) {
    auto scale = idx.getScale();
    if (scale < 0) MSG_ABORT("Negative value in bit-shift");
    int two_n = 1 << (scale + 1);

    int l[D];
    for (auto i = 0; i < D; i++) {
        l[i] = idx[i] + two_n / 2;
        if (l[i] >= two_n) return false;
        if (l[i] < 0) return false;
    }

    return true;
}

template <int D> void index_manipulation(NodeIndex<D> &idx, const std::array<bool, D> &periodic) {
    if (not periodic[0]) MSG_ABORT("Only for periodic cases!");

    auto scale = idx.getScale();
    if (scale < 0) {
        std::array<int, D> translation;
        for (auto i = 0; i < D; i++) {
            if (idx[i] < 0) translation[i] = -1;
            if (idx[i] >= 0) translation[i] = 0;
        }
        idx.setTranslation(translation);

    } else {
        std::array<int, D> translation;
        int two_n = 1 << (scale + 1);

        for (auto i = 0; i < D; i++) {
            translation[i] = idx[i] + two_n / 2;
            if (periodic[i]) {
                if (translation[i] >= two_n) translation[i] = translation[i] % two_n;
                if (translation[i] < 0) translation[i] = (translation[i] + 1) % two_n + two_n - 1;
            }
            translation[i] -= two_n / 2;
        }
        idx.setTranslation(translation);
    }
}

// Assumes a coordiate for a function defined on the [-1, 1] periodic cell.
// If r[i] is outside the unit-cell the function value is mapped to the unit-cell.
template <int D> void coord_manipulation(Coord<D> &r, const std::array<bool, D> &periodic) {
    for (auto i = 0; i < D; i++) {
        r[i] *= 0.5;
        r[i] += 0.5;
        if (periodic[i]) {
            if (r[i] >= 1.0) r[i] = std::fmod(r[i], 1.0);
            if (r[i] < 0.0) r[i] = std::fmod(r[i], 1.0) + 1.0;
        }
        r[i] -= 0.5;
        r[i] /= 0.5;
    }
}

template bool in_unit_cell<1>(NodeIndex<1> idx);
template bool in_unit_cell<2>(NodeIndex<2> idx);
template bool in_unit_cell<3>(NodeIndex<3> idx);

template void index_manipulation<1>(NodeIndex<1> &idx, const std::array<bool, 1> &periodic);
template void index_manipulation<2>(NodeIndex<2> &idx, const std::array<bool, 2> &periodic);
template void index_manipulation<3>(NodeIndex<3> &idx, const std::array<bool, 3> &periodic);

template void coord_manipulation<1>(Coord<1> &r, const std::array<bool, 1> &periodic);
template void coord_manipulation<2>(Coord<2> &r, const std::array<bool, 2> &periodic);
template void coord_manipulation<3>(Coord<3> &r, const std::array<bool, 3> &periodic);

} // namespace periodic
} // namespace mrcpp
