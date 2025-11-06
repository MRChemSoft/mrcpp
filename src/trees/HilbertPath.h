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

/**
 * @file HilbertPath.h
 * @brief Lookup-based helper to traverse octree/quadtree children in Hilbert order.
 *
 * @details
 * A Hilbert curve traversal depends on an **orientation state** that changes
 * from parent to child. This lightweight class stores the current state and
 * provides constant-time conversions between:
 * - the **Hilbert child index** \f$h \in \{0,\dots,2^D-1\}\f$ for the
 *   current state, and
 * - the corresponding **Z-order (Morton) index** \f$z\f$,
 * as well as the **next orientation state** after descending to child \f$h\f$.
 *
 * The mappings are implemented via static lookup tables (declared here, defined
 * in the corresponding translation unit). The template parameter @p D is the
 * spatial dimension; typical values are 2 (quadtree) or 3 (octree).
 */

/**
 * @class HilbertPath
 * @tparam D Spatial dimension (e.g., 2 for quadtree, 3 for octree).
 *
 * @brief Traverse the leaf nodes of a tree following the Hilbert space-filling curve.
 *
 * @details
 * Each node visit in a Hilbert traversal has an associated **state** that
 * determines how the children are ordered. Given the current state:
 * - @ref getZIndex maps a Hilbert child index to the corresponding Morton
 *   (Z-order) child index;
 * - @ref getHIndex performs the inverse mapping (Morton to Hilbert); and
 * - @ref getChildPath returns the orientation state to use after descending
 *   to a specific Hilbert child.
 */
template <int D>
class HilbertPath final {
public:
    /** 
     * @brief Default constructor 
     */
    HilbertPath() = default;
    /** 
     * @brief Copy constructor 
     */
    HilbertPath(const HilbertPath<D> &p)
            : path(p.path) {}
    /**
     * @brief Construct a child path from a parent path and a child index.
     *
     * @param[in] p    Parent @ref HilbertPath state.
     * @param[in] cIdx Child index expressed in **Morton (Z-order)** for this parent.
     *
     * @details
     * The provided @p cIdx is first converted to the corresponding **Hilbert**
     * index for the parent state, then the next orientation state is selected
     * via the transition table.
     */
    HilbertPath(const HilbertPath<D> &p, int cIdx) {
        int hIdx = p.getHIndex(cIdx);
        this->path = p.getChildPath(hIdx);
    }
    /** 
     * @brief Assignment operator 
     */
    HilbertPath &operator=(const HilbertPath<D> &p) {
        this->path = p.path;
        return *this;
    }
    short int getPath() const { return this->path; } ///< @return the current path */
    /**
     * @brief Get path index of selected child
     *
     * @param hIdx Child index in **Hilbert** order for the current state.
     * @return Path index for the selected child
     */
    short int getChildPath(int hIdx) const { return this->pTable[this->path][hIdx]; }
    /**
     * @brief Map Hilbert child index to Morton (Z-order) child index
     *
     * @param hIdx Child index in **Hilbert** order
     * @return **Morton** child index.
     */
    int getZIndex(int hIdx) const { return this->zTable[this->path][hIdx]; }

    /**
     * @brief Map Morton (Z-order) child index to Hilbert child index.
     *
     * @param zIdx Child index in **Morton** order
     * @return **Hilbert** child index
     */
    int getHIndex(int zIdx) const { return this->hTable[this->path][zIdx]; }

private:
    /// Current Hilbert orientation state (table row selector).
    short int path{0};

    /**
     * @name Lookup tables (declared in header, defined in the .cpp)
     * Each table has 2^D columns (up to 8 for D=3) and one row per state.
     * 
     */
    static const short int pTable[][8]; ///< Next-state table: state × h -> state'
    static const int zTable[][8];       ///< Mapping: state × h -> z
    static const int hTable[][8];       ///< Mapping: state × z -> h
    /** @} */
};

} // namespace mrcpp