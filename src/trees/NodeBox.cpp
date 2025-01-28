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
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include "NodeBox.h"
#include "MWNode.h"
#include "utils/Printer.h"

namespace mrcpp {

template <int D, typename T>
NodeBox<D, T>::NodeBox(const NodeIndex<D> &idx, const std::array<int, D> &nb)
        : BoundingBox<D>(idx, nb)
        , nOccupied(0)
        , nodes(nullptr) {
    allocNodePointers();
}

template <int D, typename T>
NodeBox<D, T>::NodeBox(const BoundingBox<D> &box)
        : BoundingBox<D>(box)
        , nOccupied(0)
        , nodes(nullptr) {
    allocNodePointers();
}

template <int D, typename T>
NodeBox<D, T>::NodeBox(const NodeBox<D, T> &box)
        : BoundingBox<D>(box)
        , nOccupied(0)
        , nodes(nullptr) {
    allocNodePointers();
}

template <int D, typename T> void NodeBox<D, T>::allocNodePointers() {
    assert(this->nodes == nullptr);
    int nNodes = this->size();
    this->nodes = new MWNode<D, T> *[nNodes];
    for (int n = 0; n < nNodes; n++) { this->nodes[n] = nullptr; }
    this->nOccupied = 0;
}

template <int D, typename T> NodeBox<D, T>::~NodeBox() {
    deleteNodes();
}

template <int D, typename T> void NodeBox<D, T>::deleteNodes() {
    if (this->nodes == nullptr) { return; }
    for (int n = 0; n < this->size(); n++) { clearNode(n); }
    delete[] this->nodes;
    this->nodes = nullptr;
}

template <int D, typename T> void NodeBox<D, T>::setNode(int bIdx, MWNode<D, T> **node) {
    assert(bIdx >= 0);
    assert(bIdx < this->totBoxes);
    clearNode(bIdx);
    this->nodes[bIdx] = *node;
    this->nOccupied++;
    assert(this->nOccupied > 0);
    *node = nullptr;
}

template <int D, typename T> MWNode<D, T> &NodeBox<D, T>::getNode(NodeIndex<D> nIdx) {
    int bIdx = this->getBoxIndex(nIdx);
    return getNode(bIdx);
}

template <int D, typename T> MWNode<D, T> &NodeBox<D, T>::getNode(Coord<D> r) {
    int bIdx = this->getBoxIndex(r);
    if (bIdx < 0) MSG_ERROR("Coord out of bounds");
    return getNode(bIdx);
}

template <int D, typename T> MWNode<D, T> &NodeBox<D, T>::getNode(int bIdx) {
    assert(bIdx >= 0);
    assert(bIdx < this->totBoxes);
    assert(this->nodes[bIdx] != nullptr);
    return *this->nodes[bIdx];
}

template <int D, typename T> const MWNode<D, T> &NodeBox<D, T>::getNode(NodeIndex<D> nIdx) const {
    int bIdx = this->getBoxIndex(nIdx);
    return getNode(bIdx);
}

template <int D, typename T> const MWNode<D, T> &NodeBox<D, T>::getNode(Coord<D> r) const {
    int bIdx = this->getBoxIndex(r);
    if (bIdx < 0) MSG_ERROR("Coord out of bounds");
    return getNode(bIdx);
}

template <int D, typename T> const MWNode<D, T> &NodeBox<D, T>::getNode(int bIdx) const {
    assert(bIdx >= 0);
    assert(bIdx < this->totBoxes);
    assert(this->nodes[bIdx] != nullptr);
    return *this->nodes[bIdx];
}

template class NodeBox<1, double>;
template class NodeBox<2, double>;
template class NodeBox<3, double>;

template class NodeBox<1, ComplexDouble>;
template class NodeBox<2, ComplexDouble>;
template class NodeBox<3, ComplexDouble>;

} // namespace mrcpp
