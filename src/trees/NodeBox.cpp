/*
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include "NodeBox.h"
#include "MWNode.h"
#include "utils/Printer.h"

using namespace std;

namespace mrcpp {

template<int D>
NodeBox<D>::NodeBox() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
NodeBox<D>::NodeBox(const NodeIndex<D> &idx, const int *nb)
        : BoundingBox<D>(idx, nb),
          nOccupied(0),
          nodes(0) {
    allocNodePointers();
}

template<int D>
NodeBox<D>::NodeBox(const BoundingBox<D> &box)
        : BoundingBox<D>(box),
          nOccupied(0),
          nodes(0) {
    allocNodePointers();
}

template<int D>
NodeBox<D>::NodeBox(const NodeBox<D> &box)
        : BoundingBox<D>(box),
          nOccupied(0),
          nodes(0) {
    allocNodePointers();
}

template<int D>
NodeBox<D>& NodeBox<D>::operator=(const NodeBox<D> &box) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void NodeBox<D>::allocNodePointers() {
    assert(this->nodes == 0);
    int nNodes = this->size();
    this->nodes = new MWNode<D>*[nNodes];
    for (int n = 0; n < nNodes; n++) {
        this->nodes[n] = 0;
    }
    this->nOccupied = 0;
}

template<int D>
NodeBox<D>::~NodeBox() {
    deleteNodes();
}

template<int D>
void NodeBox<D>::deleteNodes() {
    if (this->nodes == 0) {
        return;
    }
    for (int n = 0; n < this->size(); n++) {
        clearNode(n);
    }
    delete [] this->nodes;
    this->nodes = 0;
}

template<int D>
void NodeBox<D>::setNode(int bIdx, MWNode<D> **node) {
    assert(bIdx >= 0);
    assert(bIdx < this->nBoxes[D]);
    clearNode(bIdx);
    this->nodes[bIdx] = *node;
    this->nOccupied++;
    assert(this->nOccupied > 0);
    *node = 0;
}

template<int D>
MWNode<D>& NodeBox<D>::getNode(const NodeIndex<D> &nIdx) {
    int bIdx = this->getBoxIndex(nIdx);
    return getNode(bIdx);
}

template<int D>
MWNode<D>& NodeBox<D>::getNode(const double *r) {
    int bIdx = this->getBoxIndex(r);
    if (bIdx < 0) MSG_ERROR("Coord out of bounds");
    return getNode(bIdx);
}

template<int D>
MWNode<D>& NodeBox<D>::getNode(int bIdx) {
    assert(bIdx >= 0);
    assert(bIdx < this->nBoxes[D]);
    assert(this->nodes[bIdx] != 0);
    return *this->nodes[bIdx];
}

template<int D>
const MWNode<D>& NodeBox<D>::getNode(const NodeIndex<D> &nIdx) const {
    int bIdx = this->getBoxIndex(nIdx);
    return getNode(bIdx);
}

template<int D>
const MWNode<D>& NodeBox<D>::getNode(const double *r) const {
    int bIdx = this->getBoxIndex(r);
    if (bIdx < 0) MSG_ERROR("Coord out of bounds");
    return getNode(bIdx);
}

template<int D>
const MWNode<D>& NodeBox<D>::getNode(int bIdx) const {
    assert(bIdx >= 0);
    assert(bIdx < this->nBoxes[D]);
    assert(this->nodes[bIdx] != 0);
    return *this->nodes[bIdx];
}

template class NodeBox<1>;
template class NodeBox<2>;
template class NodeBox<3>;

} // namespace mrcpp
