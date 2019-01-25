/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
 *
 *
 *          CTCC, University of Tromsø
 *
 */

#include "constants.h"

#include "BoundingBox.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

namespace mrcpp {

template <int D>
BoundingBox<D>::BoundingBox(int n,
                            const std::array<int, D> &l,
                            const std::array<int, D> &nb,
                            const std::array<double, D> &sf)
        : cornerIndex(n, l.data())
        , periodic(false) {
    setNBoxes(nb);
    setScalingFactor(sf);
    setDerivedParameters();
}

template <int D>
BoundingBox<D>::BoundingBox(const NodeIndex<D> &idx, const std::array<int, D> &nb, const std::array<double, D> &sf)
        : cornerIndex(idx)
        , periodic(false) {
    setNBoxes(nb);
    setScalingFactor(sf);
    setDerivedParameters();
}

template <int D>
BoundingBox<D>::BoundingBox(const std::array<double, D> &sf, bool pbc)
        : cornerIndex()
        , periodic(pbc) {
    setNBoxes();
    setScalingFactor(sf);
    setDerivedParameters();
}

template <int D>
BoundingBox<D>::BoundingBox(const BoundingBox<D> &box)
        : cornerIndex(box.cornerIndex)
        , periodic(box.periodic) {
    setNBoxes(box.nBoxes);
    setScalingFactor(box.getScalingFactor());
    setDerivedParameters();
}

template <int D> BoundingBox<D> &BoundingBox<D>::operator=(const BoundingBox<D> &box) {
    if (&box != this) {
        this->cornerIndex = box.cornerIndex;
        this->periodic = box.periodic;
        setNBoxes(box.nBoxes);
        setScalingFactor(box.getScalingFactor());
        setDerivedParameters();
    }
    return *this;
}

template <int D> void BoundingBox<D>::setNBoxes(const std::array<int, D> &nb) {
    this->totBoxes = 1;
    for (int d = 0; d < D; d++) {
        this->nBoxes[d] = (nb[d] > 0) ? nb[d] : 1;
        this->totBoxes *= this->nBoxes[d];
    }
    if (this->totBoxes > 1 and isPeriodic()) MSG_FATAL("Total number of boxes must be one for periodic worlds");
}

template <int D> void BoundingBox<D>::setDerivedParameters() {
    assert(this->totBoxes > 0);
    int scale = this->cornerIndex.getScale();
    const int *l = this->cornerIndex.getTranslation();
    for (int d = 0; d < D; d++) {
        assert(this->nBoxes[d] > 0);
        this->unitLengths[d] = this->scalingFactor[d] * std::pow(2.0, -scale);
        this->boxLengths[d] = this->unitLengths[d] * this->nBoxes[d];
        this->lowerBounds[d] = l[d] * this->unitLengths[d];
        this->upperBounds[d] = this->lowerBounds[d] + this->boxLengths[d];
    }
}

template <int D> void BoundingBox<D>::setScalingFactor(const std::array<double, D> &sf) {
    assert(this->totBoxes > 0);
    this->scalingFactor = sf;
    if (scalingFactor == std::array<double, D>{}) scalingFactor.fill(1.0);
}

// Specialized for D=1 below
template <int D> NodeIndex<D> BoundingBox<D>::getNodeIndex(int bIdx) const {
    assert(bIdx >= 0 and bIdx <= this->totBoxes);
    int l[D];
    for (int d = D - 1; d >= 0; d--) {
        int ncells = 1;
        for (int i = 0; i < d; i++) { ncells *= this->nBoxes[i]; }
        double div = bIdx / ncells;
        double iint;
        std::modf(div, &iint);
        l[d] = (int)iint;
        bIdx -= ncells * l[d];
    }

    int n = getScale();
    const int *cl = this->cornerIndex.getTranslation();
    for (int d = 0; d < D; d++) { l[d] += cl[d]; }
    NodeIndex<D> nIdx(n, l);
    return nIdx;
}

// Specialized for D=1 below
template <int D> int BoundingBox<D>::getBoxIndex(const Coord<D> &r) const {

    if (this->isPeriodic()) return 0;

    int idx[D];
    for (int d = 0; d < D; d++) {
        double x = r[d];
        if (x < this->lowerBounds[d]) return -1;
        if (x >= this->upperBounds[d]) return -1;
        double div = (x - this->lowerBounds[d]) / this->unitLengths[d];
        double iint;
        std::modf(div, &iint);
        idx[d] = (int)iint;
    }

    int bIdx = 0;
    for (int i = D - 1; i >= 0; i--) {
        int ncells = 1;
        for (int j = 0; j < i; j++) { ncells *= this->nBoxes[j]; }
        bIdx += ncells * idx[i];
    }
    return bIdx;
}

// Specialized for D=1 below
template <int D> int BoundingBox<D>::getBoxIndex(const NodeIndex<D> &nIdx) const {

    if (this->isPeriodic()) return 0;

    int n = nIdx.getScale();
    int cn = this->cornerIndex.getScale();
    const int *l = nIdx.getTranslation();
    const int *cl = this->cornerIndex.getTranslation();
    int relScale = n - cn;
    if (relScale < 0) return -1;

    int bIdx = 0;
    for (int d = D - 1; d >= 0; d--) {
        int ncells = 1;
        for (int i = 0; i < d; i++) { ncells *= this->nBoxes[i]; }
        int reqTransl = (l[d] >> relScale) - cl[d];
        if (reqTransl < 0 or reqTransl >= this->nBoxes[d]) return -1;
        bIdx += ncells * reqTransl;
    }
    assert(bIdx >= 0);
    assert(bIdx < this->size());
    return bIdx;
}

template <int D> std::ostream &BoundingBox<D>::print(std::ostream &o) const {
    int oldprec = Printer::setPrecision(5);
    o << std::fixed;
    if (isPeriodic()) { o << "                   The World is Periodic" << std::endl; }
    o << " total boxes      = " << size() << std::endl;
    o << " boxes            = [ ";
    for (int i = 0; i < D; i++) { o << std::setw(11) << size(i) << " "; }
    o << "]" << std::endl;
    o << " unit lengths     = [ ";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getUnitLength(i) << " "; }
    o << "]" << std::endl;
    o << " scaling factor   = [ ";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getScalingFactor(i) << " "; }
    o << "]" << std::endl;
    o << " lower bounds     = [ ";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getLowerBound(i) << " "; }
    o << "]" << std::endl;
    o << " upper bounds     = [ ";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getUpperBound(i) << " "; }
    o << "]" << std::endl;
    o << " total length     = [ ";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getBoxLength(i) << " "; }
    o << "]";
    o << std::scientific;
    Printer::setPrecision(oldprec);
    return o;
}

template <> int BoundingBox<1>::getBoxIndex(const Coord<1> &r) const {

    if (this->isPeriodic()) return 0;

    double x = r[0];
    if (x < this->lowerBounds[0]) return -1;
    if (x >= this->upperBounds[0]) return -1;
    double div = (x - this->lowerBounds[0]) / this->unitLengths[0];
    double iint;
    std::modf(div, &iint);
    return (int)iint;
}

template <> NodeIndex<1> BoundingBox<1>::getNodeIndex(int bIdx) const {
    int n = getScale();
    int cl = this->cornerIndex.getTranslation(0);
    int l = bIdx + cl;
    NodeIndex<1> nIdx(n, &l);
    return nIdx;
}

template <> int BoundingBox<1>::getBoxIndex(const NodeIndex<1> &nIdx) const {

    if (this->isPeriodic()) return 0;

    int n = nIdx.getScale();
    int l = nIdx.getTranslation(0);
    int cn = this->cornerIndex.getScale();
    int cl = this->cornerIndex.getTranslation(0);
    int relScale = n - cn;
    if (relScale < 0) return -1;

    int bIdx = (l >> relScale) - cl;
    if (bIdx < 0 or bIdx >= this->size()) {
        return -1;
    } else {
        return bIdx;
    }
}

template class BoundingBox<1>;
template class BoundingBox<2>;
template class BoundingBox<3>;

} // namespace mrcpp
