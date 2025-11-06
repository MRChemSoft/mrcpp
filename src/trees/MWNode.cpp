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
 *  Simple n-dimensional node
 */

#include "MWNode.h"
#include "MWTree.h"
#include "NodeAllocator.h"
#include "core/QuadratureCache.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"
#include "utils/parallel.h"
#include "utils/tree_utils.h"

using namespace Eigen;

namespace mrcpp {

template <int D, typename T>
MWNode<D, T>::MWNode()
        : tree(nullptr)
        , parent(nullptr)
        , nodeIndex()
        , hilbertPath() {
    setIsLeafNode();
    setIsLooseNode();

    clearNorms();
    for (int i = 0; i < getTDim(); i++) this->children[i] = nullptr;
    MRCPP_INIT_OMP_LOCK();
}

template <int D, typename T>
MWNode<D, T>::MWNode(MWTree<D, T> *tree, const NodeIndex<D> &idx)
        : tree(tree)
        , parent(nullptr)
        , nodeIndex(idx)
        , hilbertPath() {
    for (int i = 0; i < getTDim(); i++) this->children[i] = nullptr;
    clearNorms();
    clearIsAllocated();
    clearHasCoefs();
    MRCPP_INIT_OMP_LOCK();
}

template <int D, typename T>
MWNode<D, T>::MWNode(MWTree<D, T> *tree, int rIdx)
        : tree(tree)
        , parent(nullptr)
        , nodeIndex(tree->getRootBox().getNodeIndex(rIdx))
        , hilbertPath() {
    for (int i = 0; i < getTDim(); i++) this->children[i] = nullptr;
    clearNorms();
    clearIsAllocated();
    clearHasCoefs();
    MRCPP_INIT_OMP_LOCK();
}

template <int D, typename T>
MWNode<D, T>::MWNode(MWNode<D, T> *parent, int cIdx)
        : tree(parent->tree)
        , parent(parent)
        , nodeIndex(parent->getNodeIndex().child(cIdx))
        , hilbertPath(parent->getHilbertPath(), cIdx) {
    for (int i = 0; i < getTDim(); i++) this->children[i] = nullptr;
    clearNorms();
    clearIsAllocated();
    clearHasCoefs();
    MRCPP_INIT_OMP_LOCK();
}

template <int D, typename T>
MWNode<D, T>::MWNode(const MWNode<D, T> &node, bool allocCoef, bool SetCoef)
        : tree(node.tree)
        , parent(nullptr)
        , nodeIndex(node.nodeIndex)
        , hilbertPath(node.hilbertPath) {
    for (int i = 0; i < getTDim(); i++) this->children[i] = nullptr;
    setIsLeafNode();
    setIsLooseNode();
    if (allocCoef) {
        allocCoefs(this->getTDim(), this->getKp1_d());
        if (node.hasCoefs() and SetCoef) {
            setCoefBlock(0, node.getNCoefs(), node.getCoefs());
            if (this->getNCoefs() > node.getNCoefs()) {
                for (int i = node.getNCoefs(); i < this->getNCoefs(); i++) this->coefs[i] = 0.0;
            }
            this->setHasCoefs();
            this->calcNorms();
        } else {
            this->clearHasCoefs();
            this->clearNorms();
        }
    } else {
        clearHasCoefs();
        coefs = nullptr;
    }
    MRCPP_INIT_OMP_LOCK();
}

template <int D, typename T> MWNode<D, T>::~MWNode() {
    if (this->isLooseNode()) this->freeCoefs();
    MRCPP_DESTROY_OMP_LOCK();
}

/* Dummy deallocation of MWNode coefficients.
 *
 * This is just to make sure this method never really gets
 * called (derived classes must implement their own version). This was
 * to avoid having pure virtual methods in the base class.
 */
template <int D, typename T> void MWNode<D, T>::dealloc() {
    NOT_REACHED_ABORT;
}

template <int D, typename T> void MWNode<D, T>::allocCoefs(int n_blocks, int block_size) {
    if (this->n_coefs != 0) MSG_ABORT("n_coefs should be zero");
    if (this->isAllocated()) MSG_ABORT("Coefs already allocated");
    if (not this->isLooseNode()) MSG_ABORT("Only loose nodes here!");

    this->n_coefs = n_blocks * block_size;
    this->coefs = new T[this->n_coefs];

    this->clearHasCoefs();
    this->setIsAllocated();
}

template <int D, typename T> void MWNode<D, T>::freeCoefs() {
    if (not this->isLooseNode()) MSG_ABORT("Only loose nodes here!");

    if (this->coefs != nullptr) delete[] this->coefs;

    this->coefs = nullptr;
    this->n_coefs = 0;

    this->clearHasCoefs();
    this->clearIsAllocated();
}

template <int D, typename T> void MWNode<D, T>::printCoefs() const {
    if (not this->isAllocated()) MSG_ABORT("Node is not allocated");
    println(0, "\nMW coefs");
    int kp1_d = this->getKp1_d();
    for (int i = 0; i < this->n_coefs; i++) {
        if (i % kp1_d == 0) println(0, "\n");
        println(0, this->coefs[i]);
    }
}

template <int D, typename T> void MWNode<D, T>::getCoefs(Eigen::Matrix<T, Eigen::Dynamic, 1> &c) const {
    if (not this->isAllocated()) MSG_ABORT("Node is not allocated");
    if (not this->hasCoefs()) MSG_ABORT("Node has no coefs");
    if (this->n_coefs == 0) MSG_ABORT("ncoefs == 0");

    c = Eigen::Matrix<T, Eigen::Dynamic, 1>::Map(this->coefs, this->n_coefs);
}

template <int D, typename T> void MWNode<D, T>::zeroCoefs() {
    if (not this->isAllocated()) MSG_ABORT("Coefs not allocated " << *this);

    for (int i = 0; i < this->n_coefs; i++) { this->coefs[i] = 0.0; }
    this->zeroNorms();
    this->setHasCoefs();
}

template <int D, typename T> void MWNode<D, T>::attachCoefs(T *coefs) {
    this->coefs = coefs;
    this->setHasCoefs();
}

template <int D, typename T> void MWNode<D, T>::setCoefBlock(int block, int block_size, const T *c) {
    if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
    for (int i = 0; i < block_size; i++) { this->coefs[block * block_size + i] = c[i]; }
}

template <int D, typename T> void MWNode<D, T>::addCoefBlock(int block, int block_size, const T *c) {
    if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
    for (int i = 0; i < block_size; i++) { this->coefs[block * block_size + i] += c[i]; }
}

template <int D, typename T> void MWNode<D, T>::zeroCoefBlock(int block, int block_size) {
    if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
    for (int i = 0; i < block_size; i++) { this->coefs[block * block_size + i] = 0.0; }
}

template <int D, typename T> void MWNode<D, T>::giveChildrenCoefs(bool overwrite) {
    assert(this->isBranchNode());
    if (not this->isAllocated()) MSG_ABORT("Not allocated!");
    if (not this->hasCoefs()) MSG_ABORT("No coefficients!");

    if (overwrite) {
        for (int i = 0; i < getTDim(); i++) getMWChild(i).zeroCoefs();
    }

    // coeff of child should be have been allocated already here
    int stride = getMWChild(0).getNCoefs();
    T *inp = getCoefs();
    T *out = getMWChild(0).getCoefs();
    bool readOnlyScaling = false;
    if (this->isGenNode()) readOnlyScaling = true;

    auto &tree = getMWTree();
    tree_utils::mw_transform(tree, inp, out, readOnlyScaling, stride, overwrite);

    for (int i = 0; i < getTDim(); i++) {
        getMWChild(i).setHasCoefs();
        getMWChild(i).calcNorms(); // should need to compute only scaling norms
    }
}

template <int D, typename T> void MWNode<D, T>::giveChildCoefs(int cIdx, bool overwrite) {

    MWNode<D, T> node_i = *this;

    node_i.mwTransform(Reconstruction);

    int kp1_d = this->getKp1_d();
    int nChildren = this->getTDim();

    if (this->children[cIdx] == nullptr) MSG_ABORT("Child does not exist!");
    MWNode<D, T> &child = getMWChild(cIdx);
    if (overwrite) {
        child.setCoefBlock(0, kp1_d, &node_i.getCoefs()[cIdx * kp1_d]);
    } else {
        child.addCoefBlock(0, kp1_d, &node_i.getCoefs()[cIdx * kp1_d]);
    }
    child.setHasCoefs();
    child.calcNorms();
}

template <int D, typename T> void MWNode<D, T>::giveParentCoefs(bool overwrite) {
    MWNode<D, T> node = *this;
    MWNode<D, T> &parent = getMWParent();
    int kp1_d = this->getKp1_d();
    if (node.getScale() == 0) {
        NodeBox<D, T> &box = this->getMWTree().getRootBox();
        auto reverse = getTDim() - 1;
        for (auto i = 0; i < getTDim(); i++) { parent.setCoefBlock(i, kp1_d, &box.getNode(reverse - i).getCoefs()[0]); }
    } else {
        for (auto i = 0; i < getTDim(); i++) { parent.setCoefBlock(i, kp1_d, &node.getCoefs()[0]); }
    }
    parent.mwTransform(Compression);
    parent.setHasCoefs();
    parent.calcNorms();
}

template <int D, typename T> void MWNode<D, T>::copyCoefsFromChildren() {
    int kp1_d = this->getKp1_d();
    int nChildren = this->getTDim();
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        MWNode<D, T> &child = getMWChild(cIdx);
        if (not child.hasCoefs()) MSG_ABORT("Child has no coefs");
        setCoefBlock(cIdx, kp1_d, child.getCoefs());
    }
}

template <int D, typename T> void MWNode<D, T>::threadSafeGenChildren() {
    if (tree->isLocal) { NOT_IMPLEMENTED_ABORT; }
    MRCPP_SET_OMP_LOCK();
    if (isLeafNode()) {
        genChildren();
        giveChildrenCoefs();
    }
    MRCPP_UNSET_OMP_LOCK();
}

template <int D, typename T> void MWNode<D, T>::threadSafeCreateChildren() {
    if (tree->isLocal) { NOT_IMPLEMENTED_ABORT; }
    MRCPP_SET_OMP_LOCK();
    if (isLeafNode()) {
        createChildren(true);
        giveChildrenCoefs();
    }
    MRCPP_UNSET_OMP_LOCK();
}

template <int D, typename T> void MWNode<D, T>::cvTransform(int operation, bool firstchild) {
    int kp1 = this->getKp1();
    int kp1_dm1 = math_utils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    int nCoefs = this->getTDim() * kp1_d;

    auto sb = this->getMWTree().getMRA().getScalingBasis();
    const MatrixXd &S = sb.getCVMap(operation);
    T o_vec[nCoefs];
    T *out_vec = o_vec;
    T *in_vec = this->coefs;

    int nChildren = this->getTDim();
    if (firstchild) nChildren = 1;
    for (int i = 0; i < D; i++) {
        for (int t = 0; t < nChildren; t++) {
            T *out = out_vec + t * kp1_d;
            T *in = in_vec + t * kp1_d;
            math_utils::apply_filter(out, in, S, kp1, kp1_dm1, 0.0);
        }
        T *tmp = in_vec;
        in_vec = out_vec;
        out_vec = tmp;
    }

    const auto scaling_factor = this->getMWTree().getMRA().getWorldBox().getScalingFactors();
    double sf_prod = 1.0;
    for (const auto &s : scaling_factor) sf_prod *= s;
    if (sf_prod <= MachineZero) sf_prod = 1.0; // When there is no scaling factor

    int np1 = getScale() + 1; // we're working on scaling coefs on next scale
    double two_fac = std::pow(2.0, D * np1) / sf_prod;
    if (operation == Backward) {
        two_fac = std::sqrt(1.0 / two_fac);
    } else {
        two_fac = std::sqrt(two_fac);
    }
    if (IS_ODD(D)) {
        for (int i = 0; i < nCoefs; i++) { this->coefs[i] = two_fac * in_vec[i]; }
    } else {
        for (int i = 0; i < nCoefs; i++) { this->coefs[i] *= two_fac; }
    }
}
/* Old interpolating version, somewhat faster
template<int D, typename T>
void MWNode<D, T>::cvTransform(int operation) {
    const ScalingBasis &sf = this->getMWTree().getMRA().getScalingBasis();
    if (sf.getScalingType() != Interpol) {
        NOT_IMPLEMENTED_ABORT;
    }

    int quadratureOrder = sf.getQuadratureOrder();
    getQuadratureCache(qc);

    double two_scale = std::pow(2.0, this->getScale() + 1);
    VectorXd modWeights = qc.getWeights(quadratureOrder);
    if (operation == Forward) {
        modWeights = modWeights.array().inverse();
        modWeights *= two_scale;
        modWeights = modWeights.array().sqrt();
    } else if (operation == Backward) {
        modWeights *= 1.0/two_scale;
        modWeights = modWeights.array().sqrt();
    } else {
        MSG_ABORT("Invalid operation");
    }

    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    int kp1_p[D];
    for (int d = 0; d < D; d++) {
        kp1_p[d] = math_utils::ipow(kp1, d);
    }

    for (int m = 0; m < this->getTDim(); m++) {
        for (int p = 0; p < D; p++) {
            int n = 0;
            for (int i = 0; i < kp1_p[D - p - 1]; i++) {
                for (int j = 0; j < kp1; j++) {
                    for (int k = 0; k < kp1_p[p]; k++) {
                        this->coefs[m * kp1_d + n] *= modWeights[j];
                        n++;
                    }
                }
            }
        }
    }
}
*/

template <int D, typename T> void MWNode<D, T>::mwTransform(int operation) {
    int kp1 = this->getKp1();
    int kp1_dm1 = math_utils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    int nCoefs = this->getTDim() * kp1_d;
    const MWFilter &filter = getMWTree().getMRA().getFilter();
    double overwrite = 0.0;

    T o_vec[nCoefs];
    T *out_vec = o_vec;
    T *in_vec = this->coefs;

    for (int i = 0; i < D; i++) {
        int mask = 1 << i;
        for (int gt = 0; gt < this->getTDim(); gt++) {
            T *out = out_vec + gt * kp1_d;
            for (int ft = 0; ft < this->getTDim(); ft++) {
                /* Operate in direction i only if the bits along other
                 * directions are identical. The bit of the direction we
                 * operate on determines the appropriate filter/operator */
                if ((gt | mask) == (ft | mask)) {
                    T *in = in_vec + ft * kp1_d;
                    int fIdx = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                    const MatrixXd &oper = filter.getSubFilter(fIdx, operation);
                    math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                    overwrite = 1.0;
                }
            }
            overwrite = 0.0;
        }
        T *tmp = in_vec;
        in_vec = out_vec;
        out_vec = tmp;
    }
    if (IS_ODD(D)) {
        for (int i = 0; i < nCoefs; i++) { this->coefs[i] = in_vec[i]; }
    }
}

template <int D, typename T> void MWNode<D, T>::clearNorms() {
    this->squareNorm = -1.0;
    for (int i = 0; i < this->getTDim(); i++) { this->componentNorms[i] = -1.0; }
}

template <int D, typename T> void MWNode<D, T>::zeroNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) { this->componentNorms[i] = 0.0; }
}

template <int D, typename T> void MWNode<D, T>::calcNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) {
        double norm_i = calcComponentNorm(i);
        this->componentNorms[i] = norm_i;
        this->squareNorm += norm_i * norm_i;
    }
}

template <int D, typename T> double MWNode<D, T>::getScalingNorm() const {
    double sNorm = this->getComponentNorm(0);
    if (sNorm >= 0.0) {
        return sNorm * sNorm;
    } else {
        return -1.0;
    }
}

template <int D, typename T> double MWNode<D, T>::getWaveletNorm() const {
    double wNorm = 0.0;
    for (int i = 1; i < this->getTDim(); i++) {
        double norm_i = this->getComponentNorm(i);
        if (norm_i >= 0.0) {
            wNorm += norm_i * norm_i;
        } else {
            wNorm = -1.0;
        }
    }
    return wNorm;
}

template <int D, typename T> double MWNode<D, T>::calcComponentNorm(int i) const {
    if (this->isGenNode() and i != 0) return 0.0;
    assert(this->isAllocated());
    assert(this->hasCoefs());

    const T *c = this->getCoefs();
    int size = this->getKp1_d();
    int start = i * size;

    double sq_norm = 0.0;
    for (int i = start; i < start + size; i++) { sq_norm += std::norm(c[i]); }
    return std::sqrt(sq_norm);
}

template <int D, typename T> void MWNode<D, T>::reCompress() {
    if (this->isGenNode()) NOT_IMPLEMENTED_ABORT;
    if (this->isBranchNode()) {
        if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
        copyCoefsFromChildren();
        mwTransform(Compression);
        this->setHasCoefs();
        this->calcNorms();
    }
}

template <int D, typename T> bool MWNode<D, T>::crop(double prec, double splitFac, bool absPrec) {
    if (this->isEndNode()) {
        return true;
    } else {
        for (int i = 0; i < this->getTDim(); i++) {
            MWNode<D, T> &child = *this->children[i];
            if (child.crop(prec, splitFac, absPrec)) {
                if (tree_utils::split_check(*this, prec, splitFac, absPrec) == false) {
                    this->deleteChildren();
                    return true;
                }
            }
        }
    }
    return false;
}

template <int D, typename T> void MWNode<D, T>::createChildren(bool coefs) {
    NOT_REACHED_ABORT;
}

template <int D, typename T> void MWNode<D, T>::genChildren() {
    NOT_REACHED_ABORT;
}

template <int D, typename T> void MWNode<D, T>::genParent() {
    NOT_REACHED_ABORT;
}

template <int D, typename T> void MWNode<D, T>::deleteChildren() {
    if (this->isLeafNode()) return;
    for (int cIdx = 0; cIdx < getTDim(); cIdx++) {
        if (this->children[cIdx] != nullptr) {
            MWNode<D, T> &child = getMWChild(cIdx);
            child.deleteChildren();
            child.dealloc();
        }
        this->children[cIdx] = nullptr;
    }
    this->childSerialIx = -1;
    this->setIsLeafNode();
}

template <int D, typename T> void MWNode<D, T>::deleteParent() {
    if (this->parent == nullptr) return;
    MWNode<D, T> &parent = getMWParent();
    parent.deleteParent();
    parent.dealloc();
    this->parentSerialIx = -1;
    this->parent = nullptr;
}

template <int D, typename T> void MWNode<D, T>::deleteGenerated() {
    if (this->isBranchNode()) {
        if (this->isEndNode()) {
            this->deleteChildren();
        } else {
            for (int cIdx = 0; cIdx < getTDim(); cIdx++) { this->getMWChild(cIdx).deleteGenerated(); }
        }
    }
}

template <int D, typename T> Coord<D> MWNode<D, T>::getCenter() const {
    auto two_n = std::pow(2.0, -getScale());
    auto scaling_factor = getMWTree().getMRA().getWorldBox().getScalingFactors();
    auto &l = getNodeIndex();
    auto r = Coord<D>{};
    for (int d = 0; d < D; d++) r[d] = scaling_factor[d] * two_n * (l[d] + 0.5);
    return r;
}

template <int D, typename T> Coord<D> MWNode<D, T>::getUpperBounds() const {
    auto two_n = std::pow(2.0, -getScale());
    auto scaling_factor = getMWTree().getMRA().getWorldBox().getScalingFactors();
    auto &l = getNodeIndex();
    auto ub = Coord<D>{};
    for (int i = 0; i < D; i++) ub[i] = scaling_factor[i] * two_n * (l[i] + 1);
    return ub;
}

template <int D, typename T> Coord<D> MWNode<D, T>::getLowerBounds() const {
    auto two_n = std::pow(2.0, -getScale());
    auto scaling_factor = getMWTree().getMRA().getWorldBox().getScalingFactors();
    auto &l = getNodeIndex();
    auto lb = Coord<D>{};
    for (int i = 0; i < D; i++) lb[i] = scaling_factor[i] * two_n * l[i];
    return lb;
}

template <int D, typename T> int MWNode<D, T>::getChildIndex(const NodeIndex<D> &nIdx) const {
    assert(isAncestor(nIdx));
    int cIdx = 0;
    int diffScale = nIdx.getScale() - getScale() - 1;
    assert(diffScale >= 0);
    for (int d = 0; d < D; d++) {
        int bit = (nIdx[d] >> (diffScale)) & 1;
        cIdx = cIdx + (bit << d);
    }
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    return cIdx;
}

template <int D, typename T> int MWNode<D, T>::getChildIndex(const Coord<D> &r) const {
    assert(hasCoord(r));
    int cIdx = 0;
    double sFac = std::pow(2.0, -getScale());
    const NodeIndex<D> &l = getNodeIndex();
    for (int d = 0; d < D; d++) {
        if (r[d] > sFac * (l[d] + 0.5)) cIdx = cIdx + (1 << d);
    }
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    return cIdx;
}

template <int D, typename T> void MWNode<D, T>::getPrimitiveQuadPts(MatrixXd &pts) const {
    int kp1 = this->getKp1();
    pts = MatrixXd::Zero(D, kp1);

    getQuadratureCache(qc);
    const VectorXd &roots = qc.getRoots(kp1);

    double sFac = std::pow(2.0, -getScale());
    const NodeIndex<D> &l = getNodeIndex();
    for (int d = 0; d < D; d++) pts.row(d) = sFac * (roots.array() + static_cast<double>(l[d]));
}

template <int D, typename T> void MWNode<D, T>::getPrimitiveChildPts(MatrixXd &pts) const {
    int kp1 = this->getKp1();
    pts = MatrixXd::Zero(D, 2 * kp1);

    getQuadratureCache(qc);
    const VectorXd &roots = qc.getRoots(kp1);

    double sFac = std::pow(2.0, -(getScale() + 1));
    const NodeIndex<D> &l = getNodeIndex();
    for (int d = 0; d < D; d++) {
        pts.row(d).segment(0, kp1) = sFac * (roots.array() + 2.0 * static_cast<double>(l[d]));
        pts.row(d).segment(kp1, kp1) = sFac * (roots.array() + 2.0 * static_cast<double>(l[d]) + 1.0);
    }
}

template <int D, typename T> void MWNode<D, T>::getExpandedQuadPts(Eigen::MatrixXd &pts) const {
    MatrixXd prim_pts;
    getPrimitiveQuadPts(prim_pts);

    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    pts = MatrixXd::Zero(D, kp1_d);

    if (D == 1) pts = prim_pts;
    if (D == 2) math_utils::tensor_expand_coords_2D(kp1, prim_pts, pts);
    if (D == 3) math_utils::tensor_expand_coords_3D(kp1, prim_pts, pts);
    if (D >= 4) NOT_IMPLEMENTED_ABORT;
}

template <int D, typename T> void MWNode<D, T>::getExpandedChildPts(MatrixXd &pts) const {
    MatrixXd prim_pts;
    getPrimitiveChildPts(prim_pts);

    int tDim = this->getTDim();
    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    pts = MatrixXd::Zero(D, tDim * kp1_d);
    MatrixXd prim_t = MatrixXd::Zero(D, kp1);
    MatrixXd exp_t = MatrixXd::Zero(D, kp1_d);

    for (int t = 0; t < tDim; t++) {
        for (int d = 0; d < D; d++) {
            int idx = (t >> d) & 1;
            prim_t.row(d) = prim_pts.block(d, idx * kp1, 1, kp1);
        }
        if (D == 1) exp_t = prim_t;
        if (D == 2) math_utils::tensor_expand_coords_2D(kp1, prim_t, exp_t);
        if (D == 3) math_utils::tensor_expand_coords_3D(kp1, prim_t, exp_t);
        if (D >= 4) NOT_IMPLEMENTED_ABORT;
        pts.block(0, t * kp1_d, D, kp1_d) = exp_t;
    }
}

template <int D, typename T> const MWNode<D, T> *MWNode<D, T>::retrieveNodeNoGen(const NodeIndex<D> &idx) const {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(this->isAncestor(idx));
    if (this->isEndNode()) { // don't return GenNodes
        return nullptr;
    }
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNodeNoGen(idx);
}

template <int D, typename T> MWNode<D, T> *MWNode<D, T>::retrieveNodeNoGen(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(this->isAncestor(idx));
    if (this->isEndNode()) { // don't return GenNodes
        return nullptr;
    }
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNodeNoGen(idx);
}

template <int D, typename T> const MWNode<D, T> *MWNode<D, T>::retrieveNodeOrEndNode(const Coord<D> &r, int depth) const {
    if (getDepth() == depth or this->isEndNode()) { return this; }
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNodeOrEndNode(r, depth);
}

template <int D, typename T> MWNode<D, T> *MWNode<D, T>::retrieveNodeOrEndNode(const Coord<D> &r, int depth) {
    if (getDepth() == depth or this->isEndNode()) { return this; }
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNodeOrEndNode(r, depth);
}

template <int D, typename T> const MWNode<D, T> *MWNode<D, T>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) const {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    // We should in principle lock before read, but it makes things slower,
    // and the EndNode status does not change (normally ;)
    if (isEndNode()) { return this; }
    int cIdx = getChildIndex(idx);
    assert(children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNodeOrEndNode(idx);
}

template <int D, typename T> MWNode<D, T> *MWNode<D, T>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    // We should in principle lock before read, but it makes things slower,
    // and the EndNode status does not change (normally ;)
    if (isEndNode()) { return this; }
    int cIdx = getChildIndex(idx);
    assert(children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNodeOrEndNode(idx);
}

template <int D, typename T> MWNode<D, T> *MWNode<D, T>::retrieveNode(const Coord<D> &r, int depth) {
    if (depth < 0) MSG_ABORT("Invalid argument");

    if (getDepth() == depth) { return this; }
    assert(hasCoord(r));
    threadSafeGenChildren();
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNode(r, depth);
}

template <int D, typename T> MWNode<D, T> *MWNode<D, T>::retrieveNode(const NodeIndex<D> &idx, bool create) {
    if (getScale() == idx.getScale()) { // we're done
        if (tree->isLocal) {
            NOT_IMPLEMENTED_ABORT;
            // has to fetch coeff in Bank. NOT USED YET
            // int ncoefs = (1 << D) * this->getKp1_d();
            // coefs = new double[ncoefs]; // TODO must be cleaned at some stage
            // coefs = new double[ncoefs]; // TODO must be cleaned at some stage
            // tree->getNodeCoeff(idx, coefs);
        }
        assert(getNodeIndex() == idx);
        return this;
    }

    assert(isAncestor(idx));
    if (create) {
        threadSafeCreateChildren();
    } else {
        threadSafeGenChildren();
    }
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNode(idx, create);
}

template <int D, typename T> MWNode<D, T> *MWNode<D, T>::retrieveParent(const NodeIndex<D> &idx) {
    if (getScale() < idx.getScale()) MSG_ABORT("Scale error")
    if (getScale() == idx.getScale()) return this;
    if (this->parent == nullptr) {
        genParent();
        giveParentCoefs();
    }
    return this->parent->retrieveParent(idx);
}

template <int D, typename T> double MWNode<D, T>::getNodeNorm(const NodeIndex<D> &idx) const {
    if (this->getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return std::sqrt(this->squareNorm);
    }
    assert(isAncestor(idx));
    if (this->isEndNode()) { // we infer norm at lower scales
        return std::sqrt(this->squareNorm * std::pow(2.0, -D * (idx.getScale() - getScale())));
    }
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->getNodeNorm(idx);
}

template <int D, typename T> bool MWNode<D, T>::hasCoord(const Coord<D> &r) const {
    double sFac = std::pow(2.0, -getScale());
    const NodeIndex<D> &l = getNodeIndex();
    //    println(1, "[" << r[0] << "," << r[1] << "," << r[2] << "]");
    //    println(1, "[" << l[0] << "," << l[1] << "," << l[2] << "]");
    //    println(1, *this);
    for (int d = 0; d < D; d++) {
        if (r[d] < sFac * l[d] or r[d] > sFac * (l[d] + 1)) {
            //            println(1, "false");
            return false;
        }
    }
    //    println(1, "true");
    return true;
}

/** Testing if nodes are compatible wrt NodeIndex and Tree (order, rootScale,
 * relPrec, etc). */
template <int D, typename T> bool MWNode<D, T>::isCompatible(const MWNode<D, T> &node) {
    NOT_IMPLEMENTED_ABORT;
    //    if (nodeIndex != node.nodeIndex) {
    //        println(0, "nodeIndex mismatch" << std::endl);
    //        return false;
    //    }
    //    if (not this->tree->checkCompatible(*node.tree)) {
    //        println(0, "tree type mismatch" << std::endl);
    //        return false;
    //    }
    //    return true;
}

template <int D, typename T> bool MWNode<D, T>::isAncestor(const NodeIndex<D> &idx) const {
    int relScale = idx.getScale() - getScale();
    if (relScale < 0) return false;
    const NodeIndex<D> &l = getNodeIndex();
    for (int d = 0; d < D; d++) {
        int reqTransl = idx[d] >> relScale;
        if (l[d] != reqTransl) return false;
    }
    return true;
}

template <int D, typename T> bool MWNode<D, T>::isDecendant(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT;
}

template <int D, typename T> std::ostream &MWNode<D, T>::print(std::ostream &o) const {
    std::string flags = "       ";
    o << getNodeIndex();
    if (isRootNode()) flags[0] = 'R';
    if (isEndNode()) flags[1] = 'E';
    if (isBranchNode()) {
        flags[2] = 'B';
    } else {
        flags[2] = 'L';
    }
    if (isGenNode()) {
        flags[3] = 'G';
    } else {
        flags[3] = 'P';
    }
    if (isAllocated()) flags[4] = 'A';
    if (hasCoefs()) flags[5] = 'C';
    o << " " << flags;
    o << " Norms (sq, s, w) = (";
    o << std::setw(12) << std::setprecision(4) << getSquareNorm() << ",";
    o << std::setw(12) << std::setprecision(4) << getScalingNorm() << ",";
    o << std::setw(12) << std::setprecision(4) << getWaveletNorm() << ")";
    return o;
}

template <int D, typename T> void MWNode<D, T>::setMaxSquareNorm() {
    auto n = this->getScale();
    this->maxWSquareNorm = calcScaledWSquareNorm();
    this->maxSquareNorm = calcScaledSquareNorm();

    if (not this->isEndNode()) {
        for (int i = 0; i < this->getTDim(); i++) {
            MWNode<D, T> &child = *this->children[i];
            child.setMaxSquareNorm();
            this->maxSquareNorm = std::max(this->maxSquareNorm, child.maxSquareNorm);
            this->maxWSquareNorm = std::max(this->maxWSquareNorm, child.maxWSquareNorm);
        }
    }
}

template <int D, typename T> void MWNode<D, T>::resetMaxSquareNorm() {
    auto n = this->getScale();
    this->maxSquareNorm = -1.0;
    this->maxWSquareNorm = -1.0;
    if (not this->isEndNode()) {
        for (int i = 0; i < this->getTDim(); i++) {
            MWNode<D, T> &child = *this->children[i];
            child.resetMaxSquareNorm();
        }
    }
}

template class MWNode<1, double>;
template class MWNode<2, double>;
template class MWNode<3, double>;
template class MWNode<1, ComplexDouble>;
template class MWNode<2, ComplexDouble>;
template class MWNode<3, ComplexDouble>;

} // namespace mrcpp
