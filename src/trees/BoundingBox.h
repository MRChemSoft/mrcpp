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

#include "NodeIndex.h"
#include "utils/details.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/** @class BoundingBox
 *
 * @brief Class defining the computational domain
 *
 * @details The computational domain is made up of a collection of D-dimensional
 * boxes on a particular length scale \f$ n \f$. The size of each box is then
 * \f$ [2^{-n}]^D \f$, i.e. higher scale means smaller boxes, and the scale
 * may be negative. The number of boxes can be different in each dimension
 * \f$ [n_x, n_y, \dots] \f$, but they must all be on the same scale (size).
 * Box translations relative to the world origin _must_ be an integer
 * multiple of the given scale size \f$ 2^{-n} \f$.
 */

template <int D> class BoundingBox {
public:
    explicit BoundingBox(std::array<int, 2> box);
    explicit BoundingBox(int n = 0, const std::array<int, D> &l = {}, const std::array<int, D> &nb = {}, const std::array<double, D> &sf = {}, bool pbc = false);
    explicit BoundingBox(const NodeIndex<D> &idx, const std::array<int, D> &nb = {}, const std::array<double, D> &sf = {});
    explicit BoundingBox(const std::array<double, D> &sf, bool pbc = true);
    BoundingBox(const std::array<double, D> &sf, std::array<bool, D> pbc);
    BoundingBox(const BoundingBox<D> &box);
    BoundingBox<D> &operator=(const BoundingBox<D> &box);
    virtual ~BoundingBox() = default;

    inline bool operator==(const BoundingBox<D> &box) const;
    inline bool operator!=(const BoundingBox<D> &box) const;

    NodeIndex<D> getNodeIndex(int bIdx) const;

    int getBoxIndex(Coord<D> r) const;
    int getBoxIndex(NodeIndex<D> nIdx) const;

    int size() const { return this->totBoxes; }
    int size(int d) const { return this->nBoxes[d]; }
    int getScale() const { return this->cornerIndex.getScale(); }
    double getScalingFactor(int d) const { return this->scalingFactor[d]; }
    double getUnitLength(int d) const { return this->unitLengths[d]; }
    double getBoxLength(int d) const { return this->boxLengths[d]; }
    double getLowerBound(int d) const { return this->lowerBounds[d]; }
    double getUpperBound(int d) const { return this->upperBounds[d]; }
    bool isPeriodic() const { return details::are_any(this->periodic, true); }
    const std::array<bool, D> &getPeriodic() const { return this->periodic; }
    const Coord<D> &getUnitLengths() const { return this->unitLengths; }
    const Coord<D> &getBoxLengths() const { return this->boxLengths; }
    const Coord<D> &getLowerBounds() const { return this->lowerBounds; }
    const Coord<D> &getUpperBounds() const { return this->upperBounds; }
    const NodeIndex<D> &getCornerIndex() const { return this->cornerIndex; }
    const std::array<double, D> &getScalingFactors() const { return this->scalingFactor; }
    friend std::ostream &operator<<(std::ostream &o, const BoundingBox<D> &box) { return box.print(o); }

protected:
    // Fundamental parameters
    NodeIndex<D> cornerIndex;    ///< Index defining the lower corner of the box
    std::array<int, D> nBoxes{}; ///< Number of boxes in each dim, last entry total
    std::array<double, D> scalingFactor{};
    std::array<bool, D> periodic{}; ///< Sets which dimension has Periodic boundary conditions.

    // Derived parameters
    int totBoxes{1};
    Coord<D> unitLengths; ///< 1/2^initialScale
    Coord<D> boxLengths;  ///< Total length (unitLength times nBoxes)
    Coord<D> lowerBounds; ///< Box lower bound (not real)
    Coord<D> upperBounds; ///< Box upper bound (not real)

    void setNBoxes(const std::array<int, D> &nb = {});
    void setDerivedParameters();
    void setScalingFactors(const std::array<double, D> &sf);
    void setPeriodic(std::array<bool, D> periodic);
    void setPeriodic(bool periodic);

    std::ostream &print(std::ostream &o) const;
};

template <int D> bool BoundingBox<D>::operator==(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return false;
    for (int d = 0; d < D; d++) {
        if (this->size(d) != box.size(d)) return false;
    }
    return true;
}

template <int D> bool BoundingBox<D>::operator!=(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return true;
    for (int d = 0; d < D; d++) {
        if (this->size(d) != box.size(d)) return true;
    }
    return false;
}

} // namespace mrcpp
