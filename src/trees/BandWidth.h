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

/**
 * @file BandWidth.h
 * @brief Lightweight storage for per-depth operator bandwidths.
 *
 * @details
 * This class stores, for each tree depth, four component band widths plus
 * a cached maximum among them. A negative width denotes “unset/empty”.
 *
 * - Rows correspond to depths in \f$\{0,\dots,\text{maxDepth}\}\f$.
 * - Columns \f$0..3\f$ are per-component widths (e.g. for blocks T, C, B, A).
 * - Column \f$4\f$ caches the **maximum** width at that depth.
 *
 * The class provides convenience accessors, mutation with automatic update of
 * the per-depth maximum, and formatted printing.
 */

#pragma once

#include <Eigen/Core>
#include <iomanip>
#include <ostream>

namespace mrcpp {

/**
 * @class BandWidth
 * @brief Container for band widths over depths and components.
 */
class BandWidth final {
public:
    /**
     * @brief Construct with storage for @p depth + 1 rows.
     * @param depth Maximum depth to allocate (inclusive).
     *
     * All entries are initialized to -1 (empty).
     */
    BandWidth(int depth = 0)
            : widths(depth + 1, 5) {
        this->clear();
    }

    /**
     * @brief Copy-construct from another instance.
     */
    BandWidth(const BandWidth &bw)
            : widths(bw.widths) {}

    /**
     * @brief Copy-assign from another instance.
     */
    BandWidth &operator=(const BandWidth &bw);

    /**
     * @brief Set all widths (including cached maxima) to -1.
     */
    void clear() { this->widths.setConstant(-1); }

    /**
     * @brief Check whether the row for @p depth is effectively empty.
     * @param depth Depth to test.
     * @return True if @p depth is out of range or the cached max is < 0.
     */
    bool isEmpty(int depth) const;

    /**
     * @brief Highest valid depth index stored.
     * @return The maximum depth (rows - 1).
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
     * @param depth Depth to query.
     * @param index Component in {0,1,2,3}.
     * @return Width for (@p depth, @p index), or -1 if @p depth is out of range.
     */
    int getWidth(int depth, int index) const { return (depth > getDepth()) ? -1 : this->widths(depth, index); }

    /**
     * @brief Set component width and update the cached per-depth maximum.
     * @param depth Depth to modify (0..getDepth()).
     * @param index Component in {0,1,2,3}.
     * @param wd    Non-negative band width.
     */
    void setWidth(int depth, int index, int wd);

    /**
     * @brief Stream pretty-printer.
     */
    friend std::ostream &operator<<(std::ostream &o, const BandWidth &bw) { return bw.print(o); }

private:
    /// Matrix of widths; columns 0..3 = components, column 4 = cached max per depth.
    Eigen::MatrixXi widths;

    /// Implementation of formatted printing.
    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp