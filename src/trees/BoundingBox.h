/**
 *
 *          CTCC, University of Tromsø
 *
 */

#pragma once

#include <iomanip>
#include <array>

#include "NodeIndex.h"

namespace mrcpp {

template<int D>
class BoundingBox {
public:
    BoundingBox(int n = 0, const std::array<int, D> &l = {}, const std::array<int, D> &nb = {}, const std::array<double, D> &sf = {});
    BoundingBox(const NodeIndex<D> &idx, const std::array<int, D> &nb = {}, const std::array<double, D> &sf = {});
    BoundingBox(const std::array<double, D> &sf, bool pbc=true);
    BoundingBox(const BoundingBox<D> &box);
    BoundingBox<D> &operator=(const BoundingBox<D> &box);
    virtual ~BoundingBox() = default;

    inline bool operator==(const BoundingBox<D> &box) const;
    inline bool operator!=(const BoundingBox<D> &box) const;

    NodeIndex<D> getNodeIndex(int bIdx) const;

    int getBoxIndex(const Coord<D> &r) const;
    int getBoxIndex(const NodeIndex<D> &nIdx) const;

    int size() const { return this->totBoxes; }
    int size(int d) const { return this->nBoxes[d]; }
    int getScale() const { return this->cornerIndex.getScale(); }
    double getScalingFactor(int d) const { return this->scalingFactor[d]; }
    double getUnitLength(int d) const { return this->unitLengths[d]; }
    double getBoxLength(int d) const { return this->boxLengths[d]; }
    double getLowerBound(int d) const { return this->lowerBounds[d]; }
    double getUpperBound(int d) const { return this->upperBounds[d]; }
    bool isPeriodic() const { return this->periodic; }
    const Coord<D> &getUnitLengths() const { return this->unitLengths; }
    const Coord<D> &getBoxLengths() const { return this->boxLengths; }
    const Coord<D> &getLowerBounds() const { return this->lowerBounds; }
    const Coord<D> &getUpperBounds() const { return this->upperBounds; }
    const NodeIndex<D> &getCornerIndex() const { return this->cornerIndex; }
    const std::array<double, D> &getScalingFactor() const {return this->scalingFactor; }
    friend std::ostream& operator<<(std::ostream &o, const BoundingBox<D> &box) { return box.print(o); }

protected:
    // Fundamental parameters
    NodeIndex<D> cornerIndex;    ///< Index defining the lower corner of the box
    std::array<int, D> nBoxes{}; ///< Number of boxes in each dim, last entry total
    std::array<double, D> scalingFactor{};
    bool periodic;

    // Derived parameters
    int totBoxes{1};
    Coord<D> unitLengths;    ///< 1/2^initialScale
    Coord<D> boxLengths;     ///< Total length (unitLength times nBoxes)
    Coord<D> lowerBounds;    ///< Box lower bound (not real)
    Coord<D> upperBounds;    ///< Box upper bound (not real)

    void setNBoxes(const std::array<int, D> &nb = {});
    void setDerivedParameters();
    void setScalingFactor(const std::array<double, D> &sf);

    std::ostream& print(std::ostream &o) const;
};

template<int D>
bool BoundingBox<D>::operator==(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return false;
    for (int d = 0; d < D; d++) {
        if (this->size(d) != box.size(d)) return false;
    }
    return true;
}

template<int D>
bool BoundingBox<D>::operator!=(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return true;
    for (int d = 0; d < D; d++) {
        if (this->size(d) != box.size(d)) return true;
    }
    return false;
}

} // namespace mrcpp
