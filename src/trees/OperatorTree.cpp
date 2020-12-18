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

#include "OperatorTree.h"
#include "BandWidth.h"
#include "LebesgueIterator.h"
#include "OperatorNode.h"
#include "OperatorNodeAllocator.h"
#include "NodeAllocator.h"
#include "utils/Printer.h"
#include "utils/tree_utils.h"

using namespace Eigen;

namespace mrcpp {

OperatorTree::OperatorTree(const MultiResolutionAnalysis<2> &mra, double np)
        : MWTree<2>(mra)
        , normPrec(np)
        , bandWidth(nullptr)
        , nodePtrStore(nullptr)
        , nodePtrAccess(nullptr) {
    if (this->normPrec < 0.0) MSG_ABORT("Negative prec");

    this->nodeAllocator_p = new OperatorNodeAllocator(this);
    this->nodeAllocator_p->allocRoots(*this);
    this->resetEndNodeTable();
}

OperatorTree::~OperatorTree() {
    clearOperNodeCache();
    clearBandWidth();
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<2> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.dealloc();
        this->rootBox.clearNode(i);
    }
    delete this->nodeAllocator_p;
}

void OperatorTree::clearBandWidth() {
    if (this->bandWidth != nullptr) delete this->bandWidth;
    this->bandWidth = nullptr;
}

void OperatorTree::calcBandWidth(double prec) {
    if (this->bandWidth != nullptr) MSG_ERROR("Band width not properly cleared");
    this->bandWidth = new BandWidth(getDepth());

    VectorXi max_transl;
    getMaxTranslations(max_transl);

    if (prec < 0.0) prec = this->normPrec;
    for (int depth = 0; depth < this->getDepth(); depth++) {
        int n = getRootScale() + depth;
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

void OperatorTree::getMaxTranslations(VectorXi &maxTransl) {
    int nScales = this->nodesAtDepth.size();
    maxTransl = VectorXi::Zero(nScales);
    LebesgueIterator<2> it(this);
    while (it.next()) {
        int n = it.getNode().getDepth();
        const NodeIndex<2> &l = it.getNode().getNodeIndex();
        maxTransl[n] = std::max(maxTransl[n], std::abs(l[0]));
        maxTransl[n] = std::max(maxTransl[n], std::abs(l[1]));
    }
}

/** Make 1D lists, adressable from [-l, l] scale by scale, of operator node
 * pointers for fast operator retrieval. This method is not thread safe,
 * since it projects missing operator nodes on the fly. Hence, it must NEVER
 * be called within a parallel region, or all hell will break loose. This is
 * not really a problem, but you have been warned.
 */
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

/** Regenerate all s/d-coeffs by backtransformation, starting at the bottom and
 * thus purifying all coefficients. Option to overwrite or add up existing
 * coefficients of BranchNodes (can be used after operator application).
 * Reimplementation of MWTree::mwTransform() without OMP, as calculation
 * of OperatorNorm is done using random vectors, which is non-deterministic
 * in parallel. FunctionTrees should be fine. */
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

/** Regenerate all scaling coeffs by MW transformation of existing s/w-coeffs
 * on coarser scales, starting at the rootNodes. Option to overwrite or add up
 * existing scaling coefficients (can be used after operator application).
 * Reimplementation of MWTree::mwTransform() without OMP, as calculation
 * of OperatorNorm is done using random vectors, which is non-deterministic
 * in parallel. FunctionTrees should be fine. */
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

std::ostream &OperatorTree::print(std::ostream &o) {
    o << std::endl << "*OperatorTree: " << this->name << std::endl;
    return MWTree<2>::print(o);
}

} // namespace mrcpp
