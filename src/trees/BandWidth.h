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

#include <Eigen/Core>
#include <iomanip>
#include <ostream>

namespace mrcpp {

/**
 * @class BandWidth
 * @brief Container for band widths over depths and components
 */
class BandWidth final {
public:
    /**
     * @brief Constructor with storage for @p depth + 1 rows
     * @param depth Maximum depth to allocate (inclusive)
     * @details All entries are initialized to -1 (empty).
     */
    BandWidth(int depth = 0)
            : widths(depth + 1, 5) {
        this->clear();
    }

    /**
     * @brief Copy-constructor
     * @param bw Instance to copy from
     */
    BandWidth(const BandWidth &bw)
            : widths(bw.widths) {}

    /// @brief Copy-assign from another instance.
    BandWidth &operator=(const BandWidth &bw);

    /// @brief Reset all width entries to -1 (empty).
    void clear() { this->widths.setConstant(-1); }

    /**
     * @brief Check whether the row for @p depth is effectively empty
     * @param depth Depth to test
     * @return True if @p depth is out of range or the last values of the row is < 0
     */
    bool isEmpty(int depth) const;

    /**
     * @brief Highest valid depth index stored
     * @return The maximum depth (rows - 1)
     */
    int getDepth() const { return this->widths.rows() - 1; }

    /**
     * @brief Cached maximum width for a depth.
     * @param depth Depth to query.
     * @return Max width at @p depth, or -1 if @p depth is out of range.
     */
    int getMaxWidth(int depth) const { return (depth > getDepth()) ? -1 : this->widths(depth, 4); }

    /**
     * @brief Component width accessor.
     * @param depth Depth index
     * @param index Component index
     * @return Width for (@p depth, @p index), or -1 if @p depth is out of range.
     */
    int getWidth(int depth, int index) const { return (depth > getDepth()) ? -1 : this->widths(depth, index); }

    /**
     * @brief Set component width and update the cached per-depth maximum.
     * @param depth Depth to modify (0..getDepth()).
     * @param index Component in {0, 1, 2, 3}.
     * @param wd    Non-negative band width.
     */
    void setWidth(int depth, int index, int wd);

    /// @brief Formatted printing of the BandWidth instance.
    friend std::ostream &operator<<(std::ostream &o, const BandWidth &bw) { return bw.print(o); }

private:
    /// Matrix of widths; columns 0..3 = components, column 4 = cached max per depth.
    Eigen::MatrixXi widths;

    /// @brief Formatted printing of the BandWidth instance.
    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp