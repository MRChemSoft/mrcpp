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

/*
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include "NodeBox.h"
#include "MRCPP/utils/Printer.h"
#include "MWNode.h"

namespace mrcpp {

template <int D>
NodeBox<D>::NodeBox(const NodeIndex<D> &idx, const std::array<int, D> &nb)
        : BoundingBox<D>(idx, nb)
        , nOccupied(0)
        , nodes(0) {
    allocNodePointers();
}

template <int D>
NodeBox<D>::NodeBox(const BoundingBox<D> &box)
        : BoundingBox<D>(box)
        , nOccupied(0)
        , nodes(0) {
    allocNodePointers();
}

template <int D>
NodeBox<D>::NodeBox(const NodeBox<D> &box)
        : BoundingBox<D>(box)
        , nOccupied(0)
        , nodes(0) {
    allocNodePointers();
}

template <int D> void NodeBox<D>::allocNodePointers() {
    assert(this->nodes == 0);
    int nNodes = this->size();
    this->nodes = new MWNode<D> *[nNodes];
    for (int n = 0; n < nNodes; n++) { this->nodes[n] = 0; }
    this->nOccupied = 0;
}

template <int D> NodeBox<D>::~NodeBox() {
    deleteNodes();
}

template <int D> void NodeBox<D>::deleteNodes() {
    if (this->nodes == 0) { return; }
    for (int n = 0; n < this->size(); n++) { clearNode(n); }
    delete[] this->nodes;
    this->nodes = 0;
}

template <int D> void NodeBox<D>::setNode(int bIdx, MWNode<D> **node) {
    assert(bIdx >= 0);
    assert(bIdx < this->totBoxes);
    clearNode(bIdx);
    this->nodes[bIdx] = *node;
    this->nOccupied++;
    assert(this->nOccupied > 0);
    *node = 0;
}

template <int D> MWNode<D> &NodeBox<D>::getNode(const NodeIndex<D> &nIdx) {
    int bIdx = this->getBoxIndex(nIdx);
    return getNode(bIdx);
}

template <int D> MWNode<D> &NodeBox<D>::getNode(const Coord<D> &r) {
    int bIdx = this->getBoxIndex(r);
    if (bIdx < 0) MSG_ERROR("Coord out of bounds");
    return getNode(bIdx);
}

template <int D> MWNode<D> &NodeBox<D>::getNode(int bIdx) {
    assert(bIdx >= 0);
    assert(bIdx < this->totBoxes);
    assert(this->nodes[bIdx] != 0);
    return *this->nodes[bIdx];
}

template <int D> const MWNode<D> &NodeBox<D>::getNode(const NodeIndex<D> &nIdx) const {
    int bIdx = this->getBoxIndex(nIdx);
    return getNode(bIdx);
}

template <int D> const MWNode<D> &NodeBox<D>::getNode(const Coord<D> &r) const {
    int bIdx = this->getBoxIndex(r);
    if (bIdx < 0) MSG_ERROR("Coord out of bounds");
    return getNode(bIdx);
}

template <int D> const MWNode<D> &NodeBox<D>::getNode(int bIdx) const {
    assert(bIdx >= 0);
    assert(bIdx < this->totBoxes);
    assert(this->nodes[bIdx] != 0);
    return *this->nodes[bIdx];
}

template class NodeBox<1>;
template class NodeBox<2>;
template class NodeBox<3>;

} // namespace mrcpp
