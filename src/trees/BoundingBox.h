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

#include <array>
#include <iomanip>
#include <ostream>

#include "NodeIndex.h"
#include "utils/details.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @file BoundingBox.h
 * @brief Declaration of the @ref BoundingBox domain descriptor.
 *
 * @details
 * The bounding box defines the computational “world” for multiresolution
 * trees. In \(D\) dimensions it is described by:
 * - a **corner index** (scale and integer translation),
 * - a **count of boxes** per dimension (all on the same scale),
 * - a **scaling factor** per dimension (physical unit lengths),
 * - and optional **periodic boundary conditions** per dimension.
 *
 * From these fundamental parameters, derived quantities such as unit
 * lengths, total box lengths and physical bounds are computed.
 */

/**
 * @class BoundingBox
 * @tparam D Spatial dimension (1, 2, or 3).
 * @brief Defines the \(D\)-dimensional computational domain (“world”).
 *
 * @details
 * The world is a Cartesian grid of equally-sized boxes at a given scale
 * \(n\). Each box has edge length \(2^{-n}\) in grid units, optionally
 * multiplied by a per-dimension scaling factor to reflect physical units.
 * The lower-back corner of the world is given by an integer translation
 * at the same scale. Periodicity can be enabled per dimension.
 */
template <int D> class BoundingBox {
public:
    /**
     * @brief Construct a non-periodic world from symmetric or half-open bounds.
     *
     * @param box Two integers \{lower, upper\}. Supported forms are
     *        \{0, L\} or \{-L, L\} with \(L>0\).
     *
     * @details
     * Chooses a root scale so that the per-dimension scaling factor
     * satisfies \(1 < \text{sfac} < 2\). The same bounds apply to all
     * dimensions. Periodicity is disabled.
     */
    explicit BoundingBox(std::array<int, 2> box);

    /**
     * @brief Fully-specified constructor (all dimensions share the same scale).
     *
     * @param n   Root scale (can be negative).
     * @param l   Integer translation (corner index) per dimension.
     * @param nb  Number of boxes per dimension (non-zero positives will be used; zeros mean 1).
     * @param sf  Scaling factor per dimension (non-positive entries are rejected).
     * @param pbc If true, all dimensions are periodic.
     *
     * @details
     * This is the most general constructor for rectangular worlds at a single
     * multiresolution scale. Periodicity is global (all-or-nothing).
     */
    explicit BoundingBox(int n = 0,
                         const std::array<int, D> &l = {},
                         const std::array<int, D> &nb = {},
                         const std::array<double, D> &sf = {},
                         bool pbc = false);

    /**
     * @brief Construct from a corner @ref NodeIndex and per-dimension sizes.
     *
     * @param idx Corner node index (scale and integer translation).
     * @param nb  Number of boxes per dimension.
     * @param sf  Scaling factor per dimension.
     *
     * @details
     * Periodicity is disabled. Useful when the corner is already known
     * in multiresolution units.
     */
    explicit BoundingBox(const NodeIndex<D> &idx,
                         const std::array<int, D> &nb = {},
                         const std::array<double, D> &sf = {});

    /**
     * @brief Construct periodic (all dimensions) world from scaling factors.
     *
     * @param sf  Scaling factor per dimension.
     * @param pbc If true, enables periodicity for all dimensions (default true).
     */
    explicit BoundingBox(const std::array<double, D> &sf, bool pbc = true);

    /**
     * @brief Construct world with per-dimension periodicity flags.
     *
     * @param sf  Scaling factor per dimension.
     * @param pbc Periodicity flags per dimension.
     */
    BoundingBox(const std::array<double, D> &sf, std::array<bool, D> pbc);

    /**
     * @brief Copy constructor.
     */
    BoundingBox(const BoundingBox<D> &box);

    /**
     * @brief Copy assignment.
     */
    BoundingBox<D> &operator=(const BoundingBox<D> &box);

    /// Defaulted virtual destructor.
    virtual ~BoundingBox() = default;

    /**
     * @name Equality
     * @brief Compare corner and per-dimension box counts.
     * @{
     */
    inline bool operator==(const BoundingBox<D> &box) const;
    inline bool operator!=(const BoundingBox<D> &box) const;
    /// @}

    /**
     * @brief Convert a world-box index to a @ref NodeIndex at the root scale.
     * @param bIdx Linear index of the box within the world.
     * @return Corner node index for that box.
     *
     * @note Specializations provide efficient versions for \(D=1\).
     */
    NodeIndex<D> getNodeIndex(int bIdx) const;

    /**
     * @brief Map a physical coordinate to the enclosing world-box index.
     * @param r Physical coordinate (scaled by @ref getScalingFactors()).
     * @return Linear index of the box, or -1 if outside and non-periodic.
     */
    int getBoxIndex(Coord<D> r) const;

    /**
     * @brief Map a @ref NodeIndex to the enclosing world-box index.
     * @param nIdx Node index (possibly at a finer scale).
     * @return Linear index of the box, or -1 if outside or at coarser scale.
     */
    int getBoxIndex(NodeIndex<D> nIdx) const;

    /// @name Size and scale queries
    /// @{
    int size() const { return this->totBoxes; }                 ///< Total number of boxes.
    int size(int d) const { return this->nBoxes[d]; }           ///< Number of boxes along dimension @p d.
    int getScale() const { return this->cornerIndex.getScale(); } ///< Root scale \(n\).
    /// @}

    /// @name Geometry (per-dimension)
    /// @{
    double getScalingFactor(int d) const { return this->scalingFactor[d]; } ///< Physical scaling factor.
    double getUnitLength(int d) const { return this->unitLengths[d]; }      ///< Unit length \(= \text{sfac}\cdot 2^{-n}\).
    double getBoxLength(int d) const { return this->boxLengths[d]; }        ///< Total world length along @p d.
    double getLowerBound(int d) const { return this->lowerBounds[d]; }      ///< Physical lower bound.
    double getUpperBound(int d) const { return this->upperBounds[d]; }      ///< Physical upper bound.
    /// @}

    /// @name Periodicity
    /// @{
    bool isPeriodic() const { return details::are_any(this->periodic, true); } ///< Any dimension periodic?
    const std::array<bool, D> &getPeriodic() const { return this->periodic; }  ///< Per-dimension flags.
    /// @}

    /// @name Bulk getters
    /// @{
    const Coord<D> &getUnitLengths() const { return this->unitLengths; }
    const Coord<D> &getBoxLengths() const { return this->boxLengths; }
    const Coord<D> &getLowerBounds() const { return this->lowerBounds; }
    const Coord<D> &getUpperBounds() const { return this->upperBounds; }
    const NodeIndex<D> &getCornerIndex() const { return this->cornerIndex; }
    const std::array<double, D> &getScalingFactors() const { return this->scalingFactor; }
    /// @}

    /**
     * @brief Pretty-printer (human-readable).
     */
    friend std::ostream &operator<<(std::ostream &o, const BoundingBox<D> &box) { return box.print(o); }

protected:
    // ---------------- Fundamental parameters ----------------

    NodeIndex<D> cornerIndex;          ///< Lower-corner node (scale + integer translation).
    std::array<int, D> nBoxes{};       ///< Number of boxes per dimension.
    std::array<double, D> scalingFactor{}; ///< Physical scaling factors per dimension.
    std::array<bool, D> periodic{};    ///< Periodicity flags per dimension.

    // ---------------- Derived parameters ----------------

    int totBoxes{1};       ///< Product of @ref nBoxes.
    Coord<D> unitLengths;  ///< Per-dimension unit length (\( \text{sfac}\cdot 2^{-n} \)).
    Coord<D> boxLengths;   ///< Total world length per dimension.
    Coord<D> lowerBounds;  ///< Physical lower bounds.
    Coord<D> upperBounds;  ///< Physical upper bounds.

    /**
     * @brief Set number of boxes per dimension.
     * @param nb If an entry is zero, it is treated as one.
     */
    void setNBoxes(const std::array<int, D> &nb = {});

    /**
     * @brief Compute all derived parameters from fundamentals.
     *
     * @details
     * Uses @ref cornerIndex, @ref nBoxes and @ref scalingFactor to fill
     * unit lengths, box lengths and physical bounds.
     */
    void setDerivedParameters();

    /**
     * @brief Set scaling factors per dimension, validating positivity.
     * @param sf Per-dimension scaling factors. Empty value means all ones.
     */
    void setScalingFactors(const std::array<double, D> &sf);

    /**
     * @brief Set periodicity per dimension.
     * @param periodic Flags per dimension.
     */
    void setPeriodic(std::array<bool, D> periodic);

    /**
     * @brief Set global periodicity (all-or-nothing).
     * @param periodic If true, all dimensions are periodic.
     */
    void setPeriodic(bool periodic);

    /**
     * @brief Print a formatted summary to stream @p o.
     */
    std::ostream &print(std::ostream &o) const;
};

// ---------------- Inline comparisons ----------------

/**
 * @brief Equality: same corner index and per-dimension box counts.
 */
template <int D> bool BoundingBox<D>::operator==(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return false;
    for (int d = 0; d < D; d++) {
        if (this->size(d) != box.size(d)) return false;
    }
    return true;
}

/**
 * @brief Inequality: differs in corner index or in any per-dimension box count.
 */
template <int D> bool BoundingBox<D>::operator!=(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return true;
    for (int d = 0; d < D; d++) {
        if (this->size(d) != box.size(d)) return true;
    }
    return false;
}

} // namespace mrcpp