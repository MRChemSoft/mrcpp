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

#include "MRCPP/constants.h"

#include "BoundingBox.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"

namespace mrcpp {

/** @returns New BoundingBox object
 *
 * @param[in] box: [lower, upper] bound in all dimensions
 *
 * @details Creates a box with appropriate root scale and scaling
 * factor to fit the given bounds, which applies to _all_ dimensions.
 * Root scale is chosen such that the scaling factor becomes 1 < sfac < 2.
 *
 * Limitations: Box must be _either_ [0,L] _or_ [-L,L], with L a positive integer.
 */
template <int D>
BoundingBox<D>::BoundingBox(std::array<int, 2> box) {
    if (box[1] < 1) {
        MSG_ERROR("Invalid upper bound: " << box[1]);
        box[1] = 1;
        MSG_WARN("Setting upper bound: " << box[1]);
    }
    if (!(box[0] == 0 or box[0] == -box[1])) {
        MSG_ERROR("Invalid lower bound: " << box[0]);
        box[0] = -box[1];
        MSG_WARN("Setting lower bound: " << box[0]);
    }
    int n = 0;
    double size = 1.0 * box[1];
    while (size >= 2.0) {
        size /= 2.0;
        n--;
    }
    auto l = std::array<int, D>{};
    auto nb = std::array<int, D>{};
    auto sfac = std::array<double, D>{};
    if (box[0] == 0) {
        l.fill(0);
        nb.fill(1);
    } else {
        l.fill(-1);
        nb.fill(2);
    }
    sfac.fill(size);
    this->cornerIndex = NodeIndex<D>(n, l);
    setPeriodic(false);
    setNBoxes(nb);
    setScalingFactors(sfac);
    setDerivedParameters();
}

/** @returns New BoundingBox object
 *
 * @param[in] n: Length scale, default 0
 * @param[in] l: Corner translation, default [0, 0, ...]
 * @param[in] nb: Number of boxes, default [1, 1, ...]
 * @param[in] sf: Scaling factor, default [1.0, 1.0, ...]
 */
template <int D>
BoundingBox<D>::BoundingBox(int n, const std::array<int, D> &l, const std::array<int, D> &nb, const std::array<double, D> &sf, bool pbc)
        : cornerIndex(n, l) {
    setPeriodic(pbc);
    setNBoxes(nb);
    setScalingFactors(sf);
    setDerivedParameters();
}

template <int D>
BoundingBox<D>::BoundingBox(const NodeIndex<D> &idx, const std::array<int, D> &nb, const std::array<double, D> &sf)
        : cornerIndex(idx) {
    setPeriodic(false);
    setNBoxes(nb);
    setScalingFactors(sf);
    setDerivedParameters();
}

template <int D>
BoundingBox<D>::BoundingBox(const std::array<double, D> &sf, bool pbc)
        : cornerIndex() {
    setPeriodic(pbc);
    setNBoxes();
    setScalingFactors(sf);
    setDerivedParameters();
}

template <int D>
BoundingBox<D>::BoundingBox(const std::array<double, D> &sf, std::array<bool, D> pbc)
        : cornerIndex() {
    setPeriodic(pbc);
    setNBoxes();
    setScalingFactors(sf);
    setDerivedParameters();
}

template <int D>
BoundingBox<D>::BoundingBox(const BoundingBox<D> &box)
        : cornerIndex(box.cornerIndex) {
    setPeriodic(box.periodic);
    setNBoxes(box.nBoxes);
    setScalingFactors(box.getScalingFactors());
    setDerivedParameters();
}

template <int D> BoundingBox<D> &BoundingBox<D>::operator=(const BoundingBox<D> &box) {
    if (&box != this) {
        this->cornerIndex = box.cornerIndex;
        this->periodic = box.periodic;
        setNBoxes(box.nBoxes);
        setScalingFactors(box.getScalingFactors());
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
}

template <int D> void BoundingBox<D>::setDerivedParameters() {
    assert(this->totBoxes > 0);
    const NodeIndex<D> &cIdx = this->cornerIndex;
    for (int d = 0; d < D; d++) {
        assert(this->nBoxes[d] > 0);
        this->unitLengths[d] = this->scalingFactor[d] * std::pow(2.0, -cIdx.getScale());
        this->boxLengths[d] = this->unitLengths[d] * this->nBoxes[d];
        this->lowerBounds[d] = cIdx[d] * this->unitLengths[d];
        this->upperBounds[d] = this->lowerBounds[d] + this->boxLengths[d];
    }
}

template <int D> void BoundingBox<D>::setScalingFactors(const std::array<double, D> &sf) {
    assert(this->totBoxes > 0);
    for (auto &x : sf)
        if (x <= 0.0 and sf != std::array<double, D>{}) MSG_ABORT("Non-positive scaling factor: " << x);
    this->scalingFactor = sf;
    if (scalingFactor == std::array<double, D>{}) scalingFactor.fill(1.0);
}

template <int D> void BoundingBox<D>::setPeriodic(bool pbc) {
    this->periodic.fill(pbc);
}

template <int D> void BoundingBox<D>::setPeriodic(std::array<bool, D> pbc) {
    this->periodic = pbc;
}

// Specialized for D=1 below
template <int D> NodeIndex<D> BoundingBox<D>::getNodeIndex(int bIdx) const {
    assert(bIdx >= 0 and bIdx <= this->totBoxes);
    std::array<int, D> l;
    for (int d = D - 1; d >= 0; d--) {
        int ncells = 1;
        for (int i = 0; i < d; i++) ncells *= this->nBoxes[i];
        double div = bIdx / ncells;
        double iint;
        std::modf(div, &iint);
        l[d] = static_cast<int>(iint);
        bIdx -= ncells * l[d];
    }

    const NodeIndex<D> &cIdx = this->cornerIndex;
    for (int d = 0; d < D; d++) l[d] += cIdx[d];
    return NodeIndex<D>(getScale(), l);
}

// Specialized for D=1 below
template <int D> int BoundingBox<D>::getBoxIndex(Coord<D> r) const {
    if (this->isPeriodic()) periodic::coord_manipulation<D>(r, this->getPeriodic());

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
        for (int j = 0; j < i; j++) ncells *= this->nBoxes[j];
        bIdx += ncells * idx[i];
    }
    return bIdx;
}

// Specialized for D=1 below
template <int D> int BoundingBox<D>::getBoxIndex(NodeIndex<D> nIdx) const {
    if (this->isPeriodic()) periodic::index_manipulation<D>(nIdx, this->getPeriodic());

    int n = nIdx.getScale();
    if (n < 0 and this->isPeriodic()) n = 0;

    const NodeIndex<D> &cIdx = this->cornerIndex;
    int relScale = n - cIdx.getScale();
    if (relScale < 0) return -1;

    int bIdx = 0;
    for (int d = D - 1; d >= 0; d--) {
        int ncells = 1;
        for (int i = 0; i < d; i++) { ncells *= this->nBoxes[i]; }
        int reqTransl = (nIdx[d] >> relScale) - cIdx[d];
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
    if (isPeriodic()) o << "                   The World is Periodic" << std::endl;
    o << " total boxes           : " << size() << std::endl;
    o << " boxes                 : [";
    for (int i = 0; i < D; i++) { o << std::setw(11) << size(i) << " "; }
    o << "]" << std::endl;
    o << " unit lengths          : [";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getUnitLength(i) << " "; }
    o << "]" << std::endl;
    o << " scaling factor        : [";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getScalingFactor(i) << " "; }
    o << "]" << std::endl;
    o << " lower bounds          : [";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getLowerBound(i) << " "; }
    o << "]" << std::endl;
    o << " upper bounds          : [";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getUpperBound(i) << " "; }
    o << "]" << std::endl;
    o << " total length          : [";
    for (int i = 0; i < D; i++) { o << std::setw(11) << getBoxLength(i) << " "; }
    o << "]";
    o << std::scientific;
    Printer::setPrecision(oldprec);
    return o;
}

template <> int BoundingBox<1>::getBoxIndex(Coord<1> r) const {

    if (this->isPeriodic()) periodic::coord_manipulation<1>(r, this->getPeriodic());

    double x = r[0];
    if (x < this->lowerBounds[0]) return -1;
    if (x >= this->upperBounds[0]) return -1;
    double div = (x - this->lowerBounds[0]) / this->unitLengths[0];
    double iint;
    std::modf(div, &iint);
    return static_cast<int>(iint);
}

template <> NodeIndex<1> BoundingBox<1>::getNodeIndex(int bIdx) const {
    const NodeIndex<1> &cIdx = this->cornerIndex;
    int n = cIdx.getScale();
    int l = bIdx + cIdx[0];
    return NodeIndex<1>(n, {l});
}

template <> int BoundingBox<1>::getBoxIndex(NodeIndex<1> nIdx) const {
    if (this->isPeriodic()) periodic::index_manipulation<1>(nIdx, this->getPeriodic());

    int n = nIdx.getScale();
    if (n < 0 and this->isPeriodic()) n = 0;
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
