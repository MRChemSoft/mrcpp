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

#include "BandWidth.h"
#include "utils/Printer.h"

namespace mrcpp {

BandWidth &BandWidth::operator=(const BandWidth &bw) = default;

bool BandWidth::isEmpty(int depth) const {
    if (depth > getDepth()) { return true; }
    if (this->widths(depth, 4) < 0) { return true; }
    return false;
}

void BandWidth::setWidth(int depth, int index, int wd) {
    assert(depth >= 0 and depth < getDepth());
    assert(index >= 0 and index < 4);
    assert(wd >= 0);
    this->widths(depth, index) = wd;
    if (wd > this->widths(depth, 4)) { this->widths(depth, 4) = wd; }
}

std::ostream &BandWidth::print(std::ostream &o) const {
    o << "  *BandWidths:" << std::endl;
    o << "   n      T   C   B   A  |  max " << std::endl;
    o << " -------------------------------" << std::endl;
    for (int depth = 0; depth <= getDepth(); depth++) {
        o << std::setw(4) << depth << " | ";
        o << std::setw(4) << this->widths(depth, 0);
        o << std::setw(4) << this->widths(depth, 1);
        o << std::setw(4) << this->widths(depth, 2);
        o << std::setw(4) << this->widths(depth, 3) << "  | ";
        o << std::setw(4) << this->widths(depth, 4) << std::endl;
    }
    o << std::endl;
    return o;
}

} // namespace mrcpp
