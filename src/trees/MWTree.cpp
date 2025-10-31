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

#include "MWTree.h"

#include "MWNode.h"
#include "MultiResolutionAnalysis.h"
#include "NodeAllocator.h"
#include "NodeIndex.h"
#include "TreeIterator.h"
#include "utils/Bank.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"
#include "utils/tree_utils.h"

using namespace Eigen;

namespace mrcpp {

template <int D, typename T>
MWTree<D, T>::MWTree(const MultiResolutionAnalysis<D> &mra, const std::string &n)
        : MRA(mra)
        , order(mra.getOrder())
        , kp1_d(math_utils::ipow(mra.getOrder() + 1, D))
        , name(n)
        , squareNorm(-1.0)
        , rootBox(mra.getWorldBox()) {
    this->nodesAtDepth.push_back(0);
}

template <int D, typename T> MWTree<D, T>::~MWTree() {
    this->endNodeTable.clear();
    if (this->nodesAtDepth.size() != 1) MSG_ERROR("Nodes at depth != 1 -> " << this->nodesAtDepth.size());
    if (this->nodesAtDepth[0] != 0) MSG_ERROR("Nodes at depth 0 != 0 -> " << this->nodesAtDepth[0]);
}

template <int D, typename T> void MWTree<D, T>::deleteRootNodes() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D, T> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.dealloc();
        this->rootBox.clearNode(i);
    }
}

template <int D, typename T> void MWTree<D, T>::clear() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D, T> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.clearHasCoefs();
        root.clearNorms();
    }
    this->resetEndNodeTable();
    this->clearSquareNorm();
}

template <int D, typename T> void MWTree<D, T>::calcSquareNorm(bool deep) {
    double treeNorm = 0.0;
    for (int n = 0; n < this->getNEndNodes(); n++) {
        MWNode<D, T> &node = getEndMWNode(n);
        if (deep) node.calcNorms();
        assert(node.hasCoefs());
        treeNorm += node.getSquareNorm();
    }
    this->squareNorm = treeNorm;
}

template <int D, typename T> void MWTree<D, T>::mwTransform(int type, bool overwrite) {
    switch (type) {
        case TopDown:
            mwTransformDown(overwrite);
            break;
        case BottomUp:
            if (not overwrite) NOT_IMPLEMENTED_ABORT;
            mwTransformUp();
            break;
        default:
            MSG_ABORT("Invalid wavelet transform");
    }
}

template <int D, typename T> void MWTree<D, T>::mwTransformUp() {
    std::vector<MWNodeVector<D, T>> nodeTable;
    tree_utils::make_node_table(*this, nodeTable);
#pragma omp parallel shared(nodeTable) num_threads(mrcpp_get_num_threads())
    {
        int start = nodeTable.size() - 2;
        for (int n = start; n >= 0; n--) {
            int nNodes = nodeTable[n].size();
#pragma omp for schedule(guided)
            for (int i = 0; i < nNodes; i++) {
                MWNode<D, T> &node = *nodeTable[n][i];
                if (node.isBranchNode()) { node.reCompress(); }
            }
        }
    }
}

template <int D, typename T> void MWTree<D, T>::mwTransformDown(bool overwrite) {
    std::vector<MWNodeVector<D, T>> nodeTable;
    tree_utils::make_node_table(*this, nodeTable);
#pragma omp parallel shared(nodeTable) num_threads(mrcpp_get_num_threads())
    {
        for (int n = 0; n < nodeTable.size(); n++) {
            int n_nodes = nodeTable[n].size();
#pragma omp for schedule(guided)
            for (int i = 0; i < n_nodes; i++) {
                MWNode<D, T> &node = *nodeTable[n][i];
                if (node.isBranchNode()) {
                    if (this->getRootScale() > node.getScale()) {
                        int reverse = n_nodes - 1;
                        int cIdx = this->getRootBox().getBoxIndex(node.getNodeIndex());
                        node.giveChildCoefs(reverse - cIdx, overwrite);
                    } else {
                        node.giveChildrenCoefs(overwrite);
                    }
                }
            }
        }
    }
}

template <int D, typename T> void MWTree<D, T>::setZero() {
    TreeIterator<D, T> it(*this);
    while (it.next()) {
        MWNode<D, T> &node = it.getNode();
        node.zeroCoefs();
    }
    this->squareNorm = 0.0;
}

template <int D, typename T> void MWTree<D, T>::incrementNodeCount(int scale) {
    int depth = scale - getRootScale();
    if (depth < 0) {
        int n = this->nodesAtNegativeDepth.size();
        if (-depth > n) {
            for (int i = 0; i < -depth - n; i++) { this->nodesAtNegativeDepth.push_back(0); }
        }
        this->nodesAtNegativeDepth[-depth - 1]++;
    } else {
        int n = this->nodesAtDepth.size() - 1;
        if (depth > n) {
            for (int i = 0; i < depth - n; i++) { this->nodesAtDepth.push_back(0); }
        }
        this->nodesAtDepth[depth]++;
    }
}

template <int D, typename T> void MWTree<D, T>::decrementNodeCount(int scale) {
    int depth = scale - getRootScale();
    if (depth < 0) {
        assert(-depth - 1 < this->nodesAtNegativeDepth.size());
        this->nodesAtNegativeDepth[-depth - 1]--;
        assert(this->nodesAtNegativeDepth[-depth - 1] >= 0);
        if (this->nodesAtNegativeDepth[-depth - 1] == 0 and this->nodesAtNegativeDepth.size() > 0) this->nodesAtNegativeDepth.pop_back();
    } else {
        assert(depth < this->nodesAtDepth.size());
        this->nodesAtDepth[depth]--;
        assert(this->nodesAtDepth[depth] >= 0);
        if (this->nodesAtDepth[depth] == 0 and this->nodesAtDepth.size() > 1) this->nodesAtDepth.pop_back();
    }
}

template <int D, typename T> int MWTree<D, T>::getNNodesAtDepth(int depth) const {
    int N = 0;
    if (depth < 0) {
        if (this->nodesAtNegativeDepth.size() >= -depth) N = this->nodesAtNegativeDepth[-depth];
    } else {
        if (this->nodesAtDepth.size() > depth) N = this->nodesAtDepth[depth];
    }
    return N;
}

template <int D, typename T> int MWTree<D, T>::getSizeNodes() const {
    auto nCoefs = 1ll * getNNodes() * getTDim() * getKp1_d();
    return sizeof(T) * nCoefs / 1024;
}

template <int D, typename T> const MWNode<D, T> *MWTree<D, T>::findNode(NodeIndex<D> idx) const {
    if (getRootBox().isPeriodic()) { periodic::index_manipulation<D>(idx, getRootBox().getPeriodic()); }
    int rIdx = getRootBox().getBoxIndex(idx);
    if (rIdx < 0) return nullptr;
    const MWNode<D, T> &root = this->rootBox.getNode(rIdx);
    assert(root.isAncestor(idx));
    return root.retrieveNodeNoGen(idx);
}

template <int D, typename T> MWNode<D, T> *MWTree<D, T>::findNode(NodeIndex<D> idx) {
    if (getRootBox().isPeriodic()) { periodic::index_manipulation<D>(idx, getRootBox().getPeriodic()); }
    int rIdx = getRootBox().getBoxIndex(idx);
    if (rIdx < 0) return nullptr;
    MWNode<D, T> &root = this->rootBox.getNode(rIdx);
    assert(root.isAncestor(idx));
    return root.retrieveNodeNoGen(idx);
}

template <int D, typename T> MWNode<D, T> &MWTree<D, T>::getNode(NodeIndex<D> idx, bool create) {
    if (getRootBox().isPeriodic()) periodic::index_manipulation<D>(idx, getRootBox().getPeriodic());

    MWNode<D, T> *out = nullptr;
    MWNode<D, T> &root = getRootBox().getNode(idx);
    if (idx.getScale() < getRootScale()) {
#pragma omp critical(gen_parent)
        out = root.retrieveParent(idx);
    } else {
        out = root.retrieveNode(idx, create);
    }
    return *out;
}

template <int D, typename T> MWNode<D, T> &MWTree<D, T>::getNodeOrEndNode(NodeIndex<D> idx) {
    if (getRootBox().isPeriodic()) { periodic::index_manipulation<D>(idx, getRootBox().getPeriodic()); }
    MWNode<D, T> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return *root.retrieveNodeOrEndNode(idx);
}

template <int D, typename T> const MWNode<D, T> &MWTree<D, T>::getNodeOrEndNode(NodeIndex<D> idx) const {
    if (getRootBox().isPeriodic()) { periodic::index_manipulation<D>(idx, getRootBox().getPeriodic()); }
    const MWNode<D, T> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return *root.retrieveNodeOrEndNode(idx);
}

template <int D, typename T> MWNode<D, T> &MWTree<D, T>::getNode(Coord<D> r, int depth) {
    MWNode<D, T> &root = getRootBox().getNode(r);
    if (depth >= 0) {
        return *root.retrieveNode(r, depth);
    } else {
        return *root.retrieveNodeOrEndNode(r, depth);
    }
}

template <int D, typename T> MWNode<D, T> &MWTree<D, T>::getNodeOrEndNode(Coord<D> r, int depth) {
    if (getRootBox().isPeriodic()) { periodic::coord_manipulation<D>(r, getRootBox().getPeriodic()); }
    MWNode<D, T> &root = getRootBox().getNode(r);
    return *root.retrieveNodeOrEndNode(r, depth);
}

template <int D, typename T> const MWNode<D, T> &MWTree<D, T>::getNodeOrEndNode(Coord<D> r, int depth) const {
    if (getRootBox().isPeriodic()) { periodic::coord_manipulation<D>(r, getRootBox().getPeriodic()); }
    const MWNode<D, T> &root = getRootBox().getNode(r);
    return *root.retrieveNodeOrEndNode(r, depth);
}

template <int D, typename T> MWNodeVector<D, T> *MWTree<D, T>::copyEndNodeTable() {
    auto *nVec = new MWNodeVector<D, T>;
    for (int n = 0; n < getNEndNodes(); n++) {
        MWNode<D, T> &node = getEndMWNode(n);
        nVec->push_back(&node);
    }
    return nVec;
}

template <int D, typename T> void MWTree<D, T>::resetEndNodeTable() {
    clearEndNodeTable();
    TreeIterator<D, T> it(*this, TopDown, Hilbert);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D, T> &node = it.getNode();
        if (node.isEndNode()) { this->endNodeTable.push_back(&node); }
    }
}

template <int D, typename T> int MWTree<D, T>::countBranchNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D, typename T> int MWTree<D, T>::countLeafNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D, typename T> int MWTree<D, T>::countNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D, typename T> int MWTree<D, T>::countAllocNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D, typename T> std::ostream &MWTree<D, T>::print(std::ostream &o) const {
    o << "  square norm: " << this->squareNorm << std::endl;
    o << "  root scale: " << this->getRootScale() << std::endl;
    o << "  order: " << this->order << std::endl;
    o << "  nodes: " << this->getNNodes() << std::endl;
    o << "  endNodes: " << this->endNodeTable.size() << std::endl;
    o << "  nodes per scale: " << std::endl;
    for (int i = this->nodesAtNegativeDepth.size() - 1; i >= 0; i--) { o << "    scale=" << -(i + this->getRootScale() + 1) << "  nodes=" << this->nodesAtNegativeDepth[i] << std::endl; }
    for (int i = 0; i < this->nodesAtDepth.size(); i++) { o << "    scale=" << i + this->getRootScale() << "  nodes=" << this->nodesAtDepth[i] << std::endl; }
    return o;
}

template <int D, typename T> void MWTree<D, T>::makeMaxSquareNorms() {
    NodeBox<D, T> &rBox = this->getRootBox();
    MWNode<D, T> **roots = rBox.getNodes();
    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
        roots[rIdx]->setMaxSquareNorm();
    }
}

template <int D, typename T> int MWTree<D, T>::getIx(NodeIndex<D> nIdx) {
    if (this->isLocal == false) MSG_ERROR("getIx only implemented in local representation");
    if (NodeIndex2serialIx.count(nIdx) == 0)
        return -1;
    else
        return NodeIndex2serialIx[nIdx];
}

template <int D, typename T> void MWTree<D, T>::getNodeCoeff(NodeIndex<D> nIdx, T *data) {
    assert(this->isLocal);
    int size = (1 << D) * kp1_d;
    int id = 0;
    for (int i = 0; i < D; i++) id += std::abs(nIdx.getTranslation(i));
    this->NodesCoeff->get_data(id, size, data);
}

template class MWTree<1, double>;
template class MWTree<2, double>;
template class MWTree<3, double>;

template class MWTree<1, ComplexDouble>;
template class MWTree<2, ComplexDouble>;
template class MWTree<3, ComplexDouble>;

} // namespace mrcpp