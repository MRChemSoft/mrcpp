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

#include "periodic_utils.h"
#include "trees/NodeIndex.h"

#include <cmath>

namespace mrcpp {
namespace periodic {

template <int D> void indx_manipulation(NodeIndex<D> &idx, std::array<bool, D> periodic) {
    int translation[D];
    int two_n = 1 << idx.getScale();
    for (auto i = 0; i < D; i++) {
        if (periodic[i]) {
            translation[i] = idx.getTranslation(i);
            if (translation[i] >= two_n) translation[i] = translation[i] % two_n;
            if (translation[i] < 0) translation[i] = (translation[i] + 1) % two_n + two_n - 1;
        }
    }
    idx.setTranslation(translation);
}

template <int D> void coord_mainpulation(Coord<D> &r, std::array<bool, D> periodic) {
    for (auto i = 0; i < D; i++) {
        if (periodic[i]) {
            if (r[i] > 1.0) r[i] = std::fmod(r[i], 1.0);
            if (r[i] < 0.0) r[i] = std::fmod(r[i], 1.0) + 1.0;
        }
    }
}

template void indx_manipulation<1>(NodeIndex<1> &idx, std::array<bool, 1> periodic);
template void indx_manipulation<2>(NodeIndex<2> &idx, std::array<bool, 2> periodic);
template void indx_manipulation<3>(NodeIndex<3> &idx, std::array<bool, 3> periodic);

template void coord_mainpulation<1>(Coord<1> &r, std::array<bool, 1> periodic);
template void coord_mainpulation<2>(Coord<2> &r, std::array<bool, 2> periodic);
template void coord_mainpulation<3>(Coord<3> &r, std::array<bool, 3> periodic);

} // namespace periodic
} // namespace mrcpp
