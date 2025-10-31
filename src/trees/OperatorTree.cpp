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

#include "OperatorTree.h"
#include "BandWidth.h"
#include "NodeAllocator.h"
#include "OperatorNode.h"
#include "TreeIterator.h"
#include "utils/Printer.h"
#include "utils/tree_utils.h"

using namespace Eigen;

namespace mrcpp {

OperatorTree::OperatorTree(const MultiResolutionAnalysis<2> &mra, double np, const std::string &name)
        : MWTree<2>(mra, name)
        , normPrec(np)
        , bandWidth(nullptr)
        , nodePtrStore(nullptr)
        , nodePtrAccess(nullptr) {
    if (this->normPrec < 0.0) MSG_ABORT("Negative prec");

    int nodesPerChunk = 1024;
    int coefsPerNode = this->getTDim() * this->getKp1_d();
    this->nodeAllocator_p = std::make_unique<NodeAllocator<2>>(this, nullptr, coefsPerNode, nodesPerChunk);
    this->allocRootNodes();
    this->resetEndNodeTable();
}

void OperatorTree::allocRootNodes() {
    auto &allocator = this->getNodeAllocator();
    auto &rootbox = this->getRootBox();

    int nRoots = rootbox.size();
    int sIdx = allocator.alloc(nRoots);

    auto n_coefs = allocator.getNCoefs();
    auto *coef_p = allocator.getCoef_p(sIdx);
    auto *root_p = allocator.getNode_p(sIdx);

    MWNode<2> **roots = rootbox.getNodes();
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
        // construct into allocator memory
        new (root_p) OperatorNode(this, rIdx);
        roots[rIdx] = root_p;

        root_p->serialIx = sIdx;
        root_p->parentSerialIx = -1;
        root_p->childSerialIx = -1;

        root_p->n_coefs = n_coefs;
        root_p->coefs = coef_p;
        root_p->setIsAllocated();

        root_p->setIsRootNode();
        root_p->setIsLeafNode();
        root_p->setIsEndNode();
        root_p->clearHasCoefs();

        this->incrementNodeCount(root_p->getScale());
        sIdx++;
        root_p++;
        coef_p += n_coefs;
    }
}

OperatorTree::~OperatorTree() {
    clearOperNodeCache();
    clearBandWidth();
    this->deleteRootNodes();
}

void OperatorTree::clearBandWidth() {
    if (this->bandWidth != nullptr) delete this->bandWidth;
    this->bandWidth = nullptr;
}

void OperatorTree::calcBandWidth(double prec) {
    if (this->bandWidth == nullptr) clearBandWidth();
    this->bandWidth = new BandWidth(getDepth());

    VectorXi max_transl;
    getMaxTranslations(max_transl);

    if (prec < 0.0) prec = this->normPrec;
    for (int depth = 0; depth < this->getDepth(); depth++) {
        int l = 0;
        bool done = false;
        while (not done) {
            done = true;
            MWNode<2> &node = getNode(depth, l);
            double thrs = std::max(MachinePrec, prec / (8.0 * (1 << depth)));
            for (int k = 0; k < 4; k++) {
                if (node.getComponentNorm(k) > thrs) {
                    this->bandWidth->setWidth(depth, k, l);
                    done = false;
                }
            }
            if (++l > max_transl[depth]) break;
        }
    }
    println(100, "\nOperator BandWidth" << *this->bandWidth);
}

bool OperatorTree::isOutsideBand(int oTransl, int o_depth, int idx) {
    return abs(oTransl) > this->bandWidth->getWidth(o_depth, idx);
}

void OperatorTree::removeRoughScaleNoise(int trust_scale) {
    MWNode<2> *p_rubbish;     // possibly inexact end node
    MWNode<2> *p_counterpart; // exact branch node
    for (int n = (this->getDepth() - 2 < trust_scale) ? this->getDepth() - 2 : trust_scale; n > this->getRootScale(); n--) {
        int N = 1 << n;
        for (int m = 0; m < N; m++)
            for (int l = 0; l < N; l++) {
                p_rubbish = this->findNode(NodeIndex<2>(n, {m, l}));
                if (p_rubbish != nullptr && p_rubbish->isEndNode()) {
                    for (int m1 = 0; m1 < N; m1++)
                        for (int l1 = 0; l1 < N; l1++)
                            if ((m1 - l1 == m - l) && (p_counterpart = this->findNode(NodeIndex<2>(n, {m1, l1}))) != nullptr && p_counterpart->isBranchNode()) {
                                for (int i = 0; i < p_counterpart->n_coefs; i++) p_rubbish->coefs[i] = p_counterpart->coefs[i];
                            }
                }
            }
        this->mwTransform(BottomUp);
    }
}

void OperatorTree::getMaxTranslations(VectorXi &maxTransl) {
    int nScales = this->nodesAtDepth.size();
    maxTransl = VectorXi::Zero(nScales);
    TreeIterator<2> it(*this);
    while (it.next()) {
        int n = it.getNode().getDepth();
        const NodeIndex<2> &l = it.getNode().getNodeIndex();
        maxTransl[n] = std::max(maxTransl[n], std::abs(l[0]));
        maxTransl[n] = std::max(maxTransl[n], std::abs(l[1]));
    }
}

void OperatorTree::setupOperNodeCache() {
    int nScales = this->nodesAtDepth.size();
    int rootScale = this->getRootScale();
    this->nodePtrStore = new OperatorNode **[nScales];
    this->nodePtrAccess = new OperatorNode **[nScales];
    VectorXi max_transl;
    getMaxTranslations(max_transl);
    for (int n = 0; n < nScales; n++) {
        int scale = rootScale + n;
        int n_transl = max_transl[n];
        int n_nodes = 2 * n_transl + 1;

        auto **nodes = new OperatorNode *[n_nodes];
        int j = 0;
        for (int i = n_transl; i >= 0; i--) {
            NodeIndex<2> idx(scale, {0, i});
            // Generated OperatorNodes are still OperatorNodes
            if (auto *oNode = dynamic_cast<OperatorNode *>(&MWTree<2>::getNode(idx))) {
                nodes[j] = oNode;
                j++;
            } else {
                NOT_REACHED_ABORT;
            }
        }
        for (int i = 1; i <= n_transl; i++) {
            NodeIndex<2> idx(scale, {i, 0});
            if (auto *oNode = dynamic_cast<OperatorNode *>(&MWTree<2>::getNode(idx))) {
                nodes[j] = oNode;
                j++;
            } else {
                NOT_REACHED_ABORT;
            }
        }

        this->nodePtrStore[n] = nodes;
        this->nodePtrAccess[n] = &nodes[n_transl];
    }
    this->resetEndNodeTable();
}

void OperatorTree::clearOperNodeCache() {
    if (this->nodePtrStore != nullptr) {
        for (int i = 0; i < getDepth(); i++) { delete[] this->nodePtrStore[i]; }
        delete[] this->nodePtrStore;
        delete[] this->nodePtrAccess;
    }
}

void OperatorTree::mwTransformUp() {
    std::vector<MWNodeVector<2>> nodeTable;
    tree_utils::make_node_table(*this, nodeTable);
    int start = nodeTable.size() - 2;
    for (int n = start; n >= 0; n--) {
        int nNodes = nodeTable[n].size();
        for (int i = 0; i < nNodes; i++) {
            MWNode<2> &node = *nodeTable[n][i];
            if (node.isBranchNode()) { node.reCompress(); }
        }
    }
}

void OperatorTree::mwTransformDown(bool overwrite) {
    std::vector<MWNodeVector<2>> nodeTable;
    tree_utils::make_node_table(*this, nodeTable);
    for (auto &n : nodeTable) {
        int n_nodes = n.size();
        for (int i = 0; i < n_nodes; i++) {
            MWNode<2> &node = *n[i];
            if (node.isBranchNode()) { node.giveChildrenCoefs(overwrite); }
        }
    }
}

std::ostream &OperatorTree::print(std::ostream &o) const {
    o << std::endl << "*OperatorTree: " << this->name << std::endl;
    return MWTree<2>::print(o);
}

} // namespace mrcpp