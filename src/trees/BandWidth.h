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

/*
 * BandWidth.h
 */

#pragma once

#include <Eigen/Core>
#include <iomanip>

namespace mrcpp {

class BandWidth final {
public:
    BandWidth(int depth = 0)
            : widths(depth + 1, 5) {
        this->clear();
    }
    BandWidth(const BandWidth &bw)
            : widths(bw.widths) {}
    BandWidth &operator=(const BandWidth &bw);

    void clear() { this->widths.setConstant(-1); }

    bool isEmpty(int depth) const;
    int getDepth() const { return this->widths.rows() - 1; }
    int getMaxWidth(int depth) const { return (depth > getDepth()) ? -1 : this->widths(depth, 4); }
    int getWidth(int depth, int index) const { return (depth > getDepth()) ? -1 : this->widths(depth, index); }
    void setWidth(int depth, int index, int wd);

    friend std::ostream &operator<<(std::ostream &o, const BandWidth &bw) { return bw.print(o); }

private:
    Eigen::MatrixXi widths; /// column 5 stores max width at depth

    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp
