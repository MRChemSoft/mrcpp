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

template <int D> class HilbertPath final {
public:
    HilbertPath() = default;
    HilbertPath(const HilbertPath<D> &p)
            : path(p.path) {}
    HilbertPath(const HilbertPath<D> &p, int cIdx) {
        int hIdx = p.getHIndex(cIdx);
        this->path = p.getChildPath(hIdx);
    }
    HilbertPath &operator=(const HilbertPath<D> &p) {
        this->path = p.path;
        return *this;
    }

    short int getPath() const { return this->path; }
    short int getChildPath(int hIdx) const;

    int getZIndex(int hIdx) const;
    int getHIndex(int zIdx) const;

private:
    short int path{0};
    static const short int pTable[][8];
    static const int zTable[][8];
    static const int hTable[][8];
};

} // namespace mrcpp