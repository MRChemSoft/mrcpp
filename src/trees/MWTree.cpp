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

/**
 *  \date April 20, 2010
 *  CTCC, University of Troms
 *
 */

#include "MWTree.h"

#include "MWNode.h"
#include "TreeIterator.h"
#include "MultiResolutionAnalysis.h"
#include "NodeAllocator.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"
#include "utils/tree_utils.h"

using namespace Eigen;

namespace mrcpp {

/** MWTree constructor with SerialTree storage for nodes.
 * Creates an empty tree object. Node construction and assignment of most of
 * the parameters are done in derived classes. */
template <int D>
MWTree<D>::MWTree(const MultiResolutionAnalysis<D> &mra, const std::string &n)
        : MRA(mra)
        , order(mra.getOrder())
        , kp1_d(math_utils::ipow(mra.getOrder() + 1, D))
        , name(n)
        , squareNorm(-1.0)
        , rootBox(mra.getWorldBox()) {
    this->nodesAtDepth.push_back(0);
}

/** MWTree destructor. */
template <int D> MWTree<D>::~MWTree() {
    this->endNodeTable.clear();
    if (this->nodesAtDepth.size() != 1) MSG_ERROR("Nodes at depth != 1 -> " << this->nodesAtDepth.size());
    if (this->nodesAtDepth[0] != 0) MSG_ERROR("Nodes at depth 0 != 0 -> " << this->nodesAtDepth[0]);
}

template <int D> void MWTree<D>::deleteRootNodes() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.dealloc();
        this->rootBox.clearNode(i);
    }
}

/** @brief Remove all nodes in the tree
 *
 * @details Leaves the tree inn the same state as after construction, i.e.
 * undefined function containing only root nodes without coefficients.
 * The assigned memory (nodeChunks in NodeAllocator) is NOT released,
 * but is immediately available to the new function.
 */
template <int D> void MWTree<D>::clear() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.clearHasCoefs();
        root.clearNorms();
    }
    this->resetEndNodeTable();
    this->clearSquareNorm();
}

/** Calculate the squared norm of a function represented as a tree.
 *
 * Norm is calculated using endNodes only, but if your endNodeTable is
 * incomplete (e.g. within growTree), the missing nodes must be given in the
 * input work vector. Involves an MPI reduction operation. */
template <int D> void MWTree<D>::calcSquareNorm() {
    double treeNorm = 0.0;
    for (int n = 0; n < this->getNEndNodes(); n++) {
        const MWNode<D> &node = getEndMWNode(n);
        assert(node.hasCoefs());
        treeNorm += node.getSquareNorm();
    }
    this->squareNorm = treeNorm;
}

template <int D> void MWTree<D>::mwTransform(int type, bool overwrite) {
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

/** Regenerate all s/d-coeffs by backtransformation, starting at the bottom and
 * thus purifying all coefficients. Option to overwrite or add up existing
 * coefficients of BranchNodes (can be used after operator application). */
template <int D> void MWTree<D>::mwTransformUp() {
    std::vector<MWNodeVector<D>> nodeTable;
    tree_utils::make_node_table(*this, nodeTable);
#pragma omp parallel shared(nodeTable) num_threads(mrcpp_get_num_threads())
    {
        int start = nodeTable.size() - 2;
        for (int n = start; n >= 0; n--) {
            int nNodes = nodeTable[n].size();
#pragma omp for schedule(guided)
            for (int i = 0; i < nNodes; i++) {
                MWNode<D> &node = *nodeTable[n][i];
                if (node.isBranchNode()) { node.reCompress(); }
            }
        }
    }
}

/** Regenerate all scaling coeffs by MW transformation of existing s/w-coeffs
 * on coarser scales, starting at the rootNodes. Option to overwrite or add up
 * existing scaling coefficients (can be used after operator application). */
template <int D> void MWTree<D>::mwTransformDown(bool overwrite) {
    std::vector<MWNodeVector<D>> nodeTable;
    tree_utils::make_node_table(*this, nodeTable);
#pragma omp parallel shared(nodeTable) num_threads(mrcpp_get_num_threads())
    {
        for (int n = 0; n < nodeTable.size(); n++) {
            int n_nodes = nodeTable[n].size();
#pragma omp for schedule(guided)
            for (int i = 0; i < n_nodes; i++) {
                MWNode<D> &node = *nodeTable[n][i];
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

/** @brief Set the MW coefficients to zero, fixed grid
 * @details Keeps the node structure of the tree, even though the zero function
 * is representable at depth zero. Use cropTree to remove unnecessary nodes.*/
template <int D> void MWTree<D>::setZero() {
    TreeIterator<D> it(*this);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        node.zeroCoefs();
    }
    this->squareNorm = 0.0;
}

/** Increment node counters for non-GenNodes. This routine is not thread
 * safe, and must NEVER be called outside a critical region in parallel.
 * It's way. way too expensive to lock the tree, so don't even think
 * about it. */
template <int D> void MWTree<D>::incrementNodeCount(int scale) {
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

/** Decrement node counters for non-GenNodes. This routine is not thread
 * safe, and must NEVER be called outside a critical region in parallel.
 * It's way. way too expensive to lock the tree, so don't even think
 * about it. */
template <int D> void MWTree<D>::decrementNodeCount(int scale) {
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

/** @returns Total number of nodes in the tree, at given depth
 * @param[in] depth: Tree depth to count, negative means count _all_ nodes
 */
template <int D> int MWTree<D>::getNNodesAtDepth(int depth) const {
    int N = 0;
    if (depth < 0) {
        if (this->nodesAtNegativeDepth.size() >= -depth) N = this->nodesAtNegativeDepth[-depth];
    } else {
        if (this->nodesAtDepth.size() > depth) N = this->nodesAtDepth[depth];
    }
    return N;
}

/** @returns Size of all MW coefs in the tree, in kB */
template <int D> int MWTree<D>::getSizeNodes() const {
    auto nCoefs = 1ll * getNNodes() * getTDim() * getKp1_d();
    return sizeof(double) * nCoefs / 1024;
}

/** Find and return the node with the given NodeIndex, const version.
 *
 * Recursive routine to find and return the node with a given NodeIndex.
 * This routine returns the appropriate Node, or a NULL pointer if
 * the node does not exist, or if it is a GenNode. Recursion starts at the
 * appropriate rootNode. */
template <int D> const MWNode<D> *MWTree<D>::findNode(const NodeIndex<D> &idx) const {
    auto nIdx = periodic::index_manipulation<D>(idx, getRootBox().getPeriodic());
    int rIdx = getRootBox().getBoxIndex(nIdx);
    if (rIdx < 0) return nullptr;
    const MWNode<D> &root = this->rootBox.getNode(rIdx);
    assert(root.isAncestor(nIdx));
    return root.retrieveNodeNoGen(nIdx);
}

/** Find and return the node with the given NodeIndex.
 *
 * Recursive routine to find and return the node with a given NodeIndex.
 * This routine returns the appropriate Node, or a NULL pointer if
 * the node does not exist, or if it is a GenNode. Recursion starts at the
 * appropriate rootNode. */
template <int D> MWNode<D> *MWTree<D>::findNode(const NodeIndex<D> &idx) {
    auto nIdx = periodic::index_manipulation<D>(idx, getRootBox().getPeriodic());
    int rIdx = getRootBox().getBoxIndex(nIdx);
    if (rIdx < 0) return nullptr;
    MWNode<D> &root = this->rootBox.getNode(rIdx);
    assert(root.isAncestor(nIdx));
    return root.retrieveNodeNoGen(nIdx);
}

/** Find and return the node with the given NodeIndex.
 *
 * This routine ALWAYS returns the node you ask for, and will generate nodes
 * that does not exist. Recursion starts at the appropriate rootNode and
 * decends from this.*/
template <int D> MWNode<D> &MWTree<D>::getNode(const NodeIndex<D> &idx) {
    auto nIdx = periodic::index_manipulation<D>(idx, getRootBox().getPeriodic());

    MWNode<D> *out = nullptr;
    MWNode<D> &root = getRootBox().getNode(nIdx);
    if (idx.getScale() < getRootScale()) {
#pragma omp critical(gen_parent)
        out = root.retrieveParent(nIdx);
    } else {
        out = root.retrieveNode(nIdx);
    }
    return *out;
}

/** Find and return the node with the given NodeIndex.
 *
 * This routine returns the Node you ask for, or the EndNode on
 * the path to the requested node, and will never create or return GenNodes.
 * Recursion starts at the appropriate rootNode and decends from this. */
template <int D> MWNode<D> &MWTree<D>::getNodeOrEndNode(const NodeIndex<D> &idx) {
    auto nIdx = periodic::index_manipulation<D>(idx, getRootBox().getPeriodic());
    MWNode<D> &root = getRootBox().getNode(nIdx);
    assert(root.isAncestor(nIdx));
    return *root.retrieveNodeOrEndNode(nIdx);
}

/** Find and return the node with the given NodeIndex.
 *
 * This routine returns the Node you ask for, or the EndNode on
 * the path to the requested node, and will never create or return GenNodes.
 * Recursion starts at the appropriate rootNode and decends from this. */
template <int D> const MWNode<D> &MWTree<D>::getNodeOrEndNode(const NodeIndex<D> &idx) const {
    auto nIdx = periodic::index_manipulation<D>(idx, getRootBox().getPeriodic());
    const MWNode<D> &root = getRootBox().getNode(nIdx);
    assert(root.isAncestor(nIdx));
    return *root.retrieveNodeOrEndNode(nIdx);
}

/** Find and return the node at a given depth that contains a given coordinate.
 *
 * This routine ALWAYS returns the node you ask for, and will generate nodes
 * that does not exist. Recursion starts at the appropriate rootNode and
 * decends from this. */
template <int D> MWNode<D> &MWTree<D>::getNode(const Coord<D> &r, int depth) {
    MWNode<D> &root = getRootBox().getNode(r);
    if (depth >= 0) {
        return *root.retrieveNode(r, depth);
    } else {
        return *root.retrieveNodeOrEndNode(r, depth);
    }
}

/** Find and return the node at a given depth that contains a given coordinate.
 *
 * This routine returns the Node you ask for, or the EndNode on
 * the path to the requested node, and will never create or return GenNodes.
 * Recursion starts at the appropriate rootNode and decends from this. */
template <int D> MWNode<D> &MWTree<D>::getNodeOrEndNode(const Coord<D> &r, int depth) {
    auto R = periodic::coord_manipulation<D>(r, getRootBox().getPeriodic());
    MWNode<D> &root = getRootBox().getNode(R);
    return *root.retrieveNodeOrEndNode(R, depth);
}

/** Find and return the node at a given depth that contains a given coordinate.
 *
 * This routine returns the Node you ask for, or the EndNode on
 * the path to the requested node, and will never create or return GenNodes.
 * Recursion starts at the appropriate rootNode and decends from this. */
template <int D> const MWNode<D> &MWTree<D>::getNodeOrEndNode(const Coord<D> &r, int depth) const {
    auto R = periodic::coord_manipulation<D>(r, getRootBox().getPeriodic());
    const MWNode<D> &root = getRootBox().getNode(R);
    return *root.retrieveNodeOrEndNode(R, depth);
}

template <int D> MWNodeVector<D> *MWTree<D>::copyEndNodeTable() {
    auto *nVec = new MWNodeVector<D>;
    for (int n = 0; n < getNEndNodes(); n++) {
        MWNode<D> &node = getEndMWNode(n);
        nVec->push_back(&node);
    }
    return nVec;
}

template <int D> void MWTree<D>::resetEndNodeTable() {
    clearEndNodeTable();
    TreeIterator<D> it(*this, TopDown, Hilbert);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        if (node.isEndNode()) { this->endNodeTable.push_back(&node); }
    }
}

template <int D> int MWTree<D>::countBranchNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> int MWTree<D>::countLeafNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
    //    int nNodes = 0;
    //    TreeIterator<D> it(*this);
    //    while (it.next()) {
    //        MWNode<D> &node = it.getNode();
    //        if (node.getDepth() == depth or depth < 0) {
    //            if (node.isLeafNode()) {
    //                nNodes++;
    //            }
    //        }
    //    }
    //    return nNodes;
}

/** Traverse tree and count nodes belonging to this rank. */
template <int D> int MWTree<D>::countNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
    //    TreeIterator<D> it(*this);
    //    int count = 0;
    //    while (it.next()) {
    //        MWNode<D> &node = it.getNode();
    //        if (node.isGenNode()) {
    //            continue;
    //        }
    //        if (not node.isForeign()) {
    //            count++;
    //        }
    //    }
    //    return count;
}

/** Traverse tree and count nodes with allocated coefficients. */
template <int D> int MWTree<D>::countAllocNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
    //    TreeIterator<D> it(*this);
    //    int count = 0;
    //    while (it.next()) {
    //        MWNode<D> &node = it.getNode();
    //        if (node.isGenNode()) {
    //            continue;
    //        }
    //        if (node.hasCoefs()) {
    //            count++;
    //        }
    //    }
    //    return count;
}

template <int D> std::ostream &MWTree<D>::print(std::ostream &o) const {
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

/** set values for maxSquareNorm in all nodes  */
template <int D> void MWTree<D>::makeMaxSquareNorms() {
    NodeBox<D> &rBox = this->getRootBox();
    MWNode<D> **roots = rBox.getNodes();
    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
        // recursively set value of children and descendants
        roots[rIdx]->setMaxSquareNorm();
    }
}

template class MWTree<1>;
template class MWTree<2>;
template class MWTree<3>;

} // namespace mrcpp
