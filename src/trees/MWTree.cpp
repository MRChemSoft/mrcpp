/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
#include "HilbertIterator.h"
#include "MultiResolutionAnalysis.h"
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
MWTree<D>::MWTree(const MultiResolutionAnalysis<D> &mra)
        : nThreads(mrcpp_get_max_threads())
        , MRA(mra)
        , order(mra.getOrder())
        , kp1_d(math_utils::ipow(mra.getOrder() + 1, D))
        , name("nn")
        , nNodes(0)
        , squareNorm(-1.0)
        , rootBox(mra.getWorldBox()) {
    this->nodesAtDepth.push_back(0);
    allocNodeCounters();
    MRCPP_INIT_OMP_LOCK();
}

/** MWTree destructor. */
template <int D> MWTree<D>::~MWTree() {
    this->endNodeTable.clear();
    if (this->nNodes != 0) MSG_ERROR("Node count != 0 -> " << this->nNodes);
    if (this->nodesAtDepth.size() != 1) MSG_ERROR("Nodes at depth != 1 -> " << this->nodesAtDepth.size());
    if (this->nodesAtDepth[0] != 0) MSG_ERROR("Nodes at depth 0 != 0 -> " << this->nodesAtDepth[0]);
    deleteNodeCounters();
    MRCPP_DESTROY_OMP_LOCK();
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
                if (node.isBranchNode()) { node.giveChildrenCoefs(overwrite); }
            }
        }
    }
}

/** @brief Set the MW coefficients to zero, fixed grid
 * @details Keeps the node structure of the tree, even though the zero function
 * is representable at depth zero. Use cropTree to remove unnecessary nodes.*/
template <int D> void MWTree<D>::setZero() {
    HilbertIterator<D> it(this);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        node.zeroCoefs();
    }
    this->squareNorm = 0.0;
}

template <int D> void MWTree<D>::allocNodeCounters() {
    this->nGenNodes = new int[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) { this->nGenNodes[i] = 0; }
}

template <int D> void MWTree<D>::deleteNodeCounters() {
    delete[] this->nGenNodes;
}

/** Increment node counters for non-GenNodes. This routine is not thread
 * safe, and must NEVER be called outside a critical region in parallel.
 * It's way. way too expensive to lock the tree, so don't even think
 * about it. */
template <int D> void MWTree<D>::incrementNodeCount(int scale) {
    int depth = scale - getRootScale();
    assert(depth >= 0);
    int n = this->nodesAtDepth.size() - 1;
    if (depth > n) {
        for (int i = 0; i < depth - n; i++) { this->nodesAtDepth.push_back(0); }
    }
    this->nodesAtDepth[depth]++;
    this->nNodes++;
}

/** Decrement node counters for non-GenNodes. This routine is not thread
 * safe, and must NEVER be called outside a critical region in parallel.
 * It's way. way too expensive to lock the tree, so don't even think
 * about it. */
template <int D> void MWTree<D>::decrementNodeCount(int scale) {
    int depth = scale - getRootScale();
    assert(depth >= 0);
    assert(depth < this->nodesAtDepth.size());
    this->nodesAtDepth[depth]--;
    assert(this->nodesAtDepth[depth] >= 0);
    if (this->nodesAtDepth[depth] == 0 and this->nodesAtDepth.size() > 1) { this->nodesAtDepth.pop_back(); }
    this->nNodes--;
    assert(this->nNodes >= 0);
}

/** Update GenNode counts in a safe way. Since GenNodes are created on the
 * fly, we cannot control when to update the node counters without locking
 * the whole tree. Therefore GenNodes update thread-private counters, which
 * get merged with the correct global counters in xxxNodes[0]. This method
 * should be called outside of the parallel region for performance reasons. */
template <int D> void MWTree<D>::updateGenNodeCounts() {
    MRCPP_SET_OMP_LOCK();
    for (int i = 1; i < this->nThreads; i++) {
        this->nGenNodes[0] += this->nGenNodes[i];
        this->nGenNodes[i] = 0;
    }
    assert(this->nGenNodes[0] >= 0);
    MRCPP_UNSET_OMP_LOCK();
}

/** Adds a GenNode to the count. */
template <int D> void MWTree<D>::incrementGenNodeCount() {
    int n = mrcpp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nGenNodes[n]++;
}

/** Removes a GenNode from the count. */
template <int D> void MWTree<D>::decrementGenNodeCount() {
    int n = mrcpp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nGenNodes[n]--;
}

/** @returns Total number of nodes in the tree, at given depth
 * @param[in] depth: Tree depth to count, negative means count _all_ nodes
 */
template <int D> int MWTree<D>::getNNodes(int depth) const {
    if (depth < 0) { return this->nNodes; }
    if (depth >= this->nodesAtDepth.size()) { return 0; }
    return this->nodesAtDepth[depth];
}

/** Get GenNode count. Includes an OMP reduction operation */
template <int D> int MWTree<D>::getNGenNodes() {
    updateGenNodeCounts();
    return this->nGenNodes[0];
}

/** @returns Size of all MW coefs in the tree, in kB */
template <int D> int MWTree<D>::getSizeNodes() const {
    auto nCoefs = 1ll * getNNodes() * getTDim() * getKp1_d();
    return sizeof(double) * nCoefs / 1024;
}

/** Find and return the node with the given NodeIndex, const version.
 *
 * Recursive routine to find and return the node with a given NodeIndex.
 * This routine returns the appropriate ProjectedNode, or a NULL pointer if
 * the node does not exist, or if it is a GenNode. Recursion starts at the
 * appropriate rootNode. */
template <int D> const MWNode<D> *MWTree<D>::findNode(NodeIndex<D> idx) const {
    if (getRootBox().isPeriodic()) { periodic::indx_manipulation<D>(idx, getRootBox().getPeriodic()); }
    int rIdx = getRootBox().getBoxIndex(idx);
    if (rIdx < 0) return nullptr;
    const MWNode<D> &root = this->rootBox.getNode(rIdx);
    assert(root.isAncestor(idx));
    return root.retrieveNodeNoGen(idx);
}

/** Find and return the node with the given NodeIndex.
 *
 * Recursive routine to find and return the node with a given NodeIndex.
 * This routine returns the appropriate ProjectedNode, or a NULL pointer if
 * the node does not exist, or if it is a GenNode. Recursion starts at the
 * appropriate rootNode. */
template <int D> MWNode<D> *MWTree<D>::findNode(NodeIndex<D> idx) {
    if (getRootBox().isPeriodic()) { periodic::indx_manipulation<D>(idx, getRootBox().getPeriodic()); }
    int rIdx = getRootBox().getBoxIndex(idx);
    if (rIdx < 0) return nullptr;
    MWNode<D> &root = this->rootBox.getNode(rIdx);
    assert(root.isAncestor(idx));
    return root.retrieveNodeNoGen(idx);
}

/** Find and return the node with the given NodeIndex.
 *
 * This routine ALWAYS returns the node you ask for, and will generate nodes
 * that does not exist. Recursion starts at the appropriate rootNode and
 * decends from this.*/
template <int D> MWNode<D> &MWTree<D>::getNode(NodeIndex<D> idx) {
    if (getRootBox().isPeriodic()) { periodic::indx_manipulation<D>(idx, getRootBox().getPeriodic()); }
    MWNode<D> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return *root.retrieveNode(idx);
}

/** Find and return the node with the given NodeIndex.
 *
 * This routine returns the ProjectedNode you ask for, or the EndNode on
 * the path to the requested node, and will never create or return GenNodes.
 * Recursion starts at the appropriate rootNode and decends from this. */
template <int D> MWNode<D> &MWTree<D>::getNodeOrEndNode(NodeIndex<D> idx) {
    if (getRootBox().isPeriodic()) { periodic::indx_manipulation<D>(idx, getRootBox().getPeriodic()); }
    MWNode<D> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return *root.retrieveNodeOrEndNode(idx);
}

/** Find and return the node with the given NodeIndex.
 *
 * This routine returns the ProjectedNode you ask for, or the EndNode on
 * the path to the requested node, and will never create or return GenNodes.
 * Recursion starts at the appropriate rootNode and decends from this. */
template <int D> const MWNode<D> &MWTree<D>::getNodeOrEndNode(NodeIndex<D> idx) const {
    if (getRootBox().isPeriodic()) { periodic::indx_manipulation<D>(idx, getRootBox().getPeriodic()); }
    const MWNode<D> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return *root.retrieveNodeOrEndNode(idx);
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
 * This routine returns the ProjectedNode you ask for, or the EndNode on
 * the path to the requested node, and will never create or return GenNodes.
 * Recursion starts at the appropriate rootNode and decends from this. */
template <int D> MWNode<D> &MWTree<D>::getNodeOrEndNode(Coord<D> r, int depth) {

    if (getRootBox().isPeriodic()) { periodic::coord_manipulation<D>(r, getRootBox().getPeriodic()); }

    MWNode<D> &root = getRootBox().getNode(r);
    return *root.retrieveNodeOrEndNode(r, depth);
}

/** Find and return the node at a given depth that contains a given coordinate.
 *
 * This routine returns the ProjectedNode you ask for, or the EndNode on
 * the path to the requested node, and will never create or return GenNodes.
 * Recursion starts at the appropriate rootNode and decends from this. */
template <int D> const MWNode<D> &MWTree<D>::getNodeOrEndNode(Coord<D> r, int depth) const {

    if (getRootBox().isPeriodic()) { periodic::coord_manipulation<D>(r, getRootBox().getPeriodic()); }
    const MWNode<D> &root = getRootBox().getNode(r);
    return *root.retrieveNodeOrEndNode(r, depth);
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
    HilbertIterator<D> it(this);
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
    //    HilbertIterator<D> it(this);
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

/** Traverse tree and count nodes. */
template <int D> void MWTree<D>::RecountNodes() {
    int DepthMax = 100, slen = 0;
    MWNode<D> *stack[DepthMax * 8];
    for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) {
        stack[slen++] = this->getRootBox().getNodes()[rIdx];
    }
    this->nNodes = 0;
    while (slen) {
        this->nNodes++;
        MWNode<D> *fpos = stack[--slen];
        for (int i = 0; i < fpos->getNChildren(); i++) stack[slen++] = fpos->children[i];
    }
}
/** Traverse tree and count nodes belonging to this rank. */
template <int D> int MWTree<D>::countNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
    //    HilbertIterator<D> it(this);
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
    //    HilbertIterator<D> it(this);
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

template <int D> void MWTree<D>::deleteGenerated() {
    for (int n = 0; n < getNEndNodes(); n++) { getEndMWNode(n).deleteGenerated(); }
}

template <int D> void MWTree<D>::saveTree(const std::string &file) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> void MWTree<D>::loadTree(const std::string &file) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> std::ostream &MWTree<D>::print(std::ostream &o) {
    o << "  square norm: " << this->squareNorm << std::endl;
    o << "  root scale: " << this->getRootScale() << std::endl;
    o << "  order: " << this->order << std::endl;
    o << "  nodes: " << this->getNNodes() << std::endl;
    o << "  endNodes: " << this->endNodeTable.size() << std::endl;
    o << "  genNodes: " << this->getNGenNodes() << std::endl;
    o << "  nodes per scale: " << std::endl;
    for (int i = 0; i < this->nodesAtDepth.size(); i++) {
        o << "    scale=" << i + this->getRootScale() << "  nodes=" << this->nodesAtDepth[i] << std::endl;
    }
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
