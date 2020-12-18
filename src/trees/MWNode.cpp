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
 *  Simple n-dimensional node
 */

#include "MWNode.h"
#include "GenNode.h"
#include "MWTree.h"
#include "ProjectedNode.h"
#include "NodeAllocator.h"
#include "core/QuadratureCache.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"
#include "utils/tree_utils.h"

using namespace Eigen;

namespace mrcpp {

/** MWNode default constructor.
 *  Should be used only by NodeAllocator to obtain
 *  virtual table pointers for the derived classes. */
template <int D>
MWNode<D>::MWNode()
        : tree(nullptr)
        , parent(nullptr)
        , squareNorm(-1.0)
        , maxSquareNorm(-1.0)
        , maxWSquareNorm(-1.0)
        , coefs(nullptr)
        , n_coefs(0)
        , nodeIndex()
        , hilbertPath()
        , status(0) {
    setIsLeafNode();
    setIsLooseNode();

    clearNorms();
    for (int i = 0; i < getTDim(); i++) { this->children[i] = nullptr; }
}

/** MWNode copy constructor.
 *  Creates loose nodes and copy coefs */
template <int D>
MWNode<D>::MWNode(const MWNode<D> &node)
        : tree(node.tree)
        , parent(nullptr)
        , squareNorm(-1.0)
        , maxSquareNorm(-1.0)
        , maxWSquareNorm(-1.0)
        , coefs(nullptr)
        , n_coefs(0)
        , nodeIndex(node.nodeIndex)
        , hilbertPath(node.hilbertPath)
        , status(0) {
    setIsLeafNode();
    setIsLooseNode();

    allocCoefs(this->getTDim(), this->getKp1_d());

    if (node.hasCoefs()) {
        setCoefBlock(0, node.getNCoefs(), node.getCoefs());
        if (this->getNCoefs() > node.getNCoefs()) {
            for (int i = node.getNCoefs(); i < this->getNCoefs(); i++) { this->coefs[i] = 0.0; }
        }
        this->setHasCoefs();
        this->calcNorms();
    } else {
        this->clearHasCoefs();
        this->clearNorms();
    }
    for (int i = 0; i < getTDim(); i++) { this->children[i] = nullptr; }
}

/** MWNode destructor.
 * Recursive deallocation of a node and all its decendants */
template <int D> MWNode<D>::~MWNode() {
    if (this->isLooseNode()) this->freeCoefs();
}

template <int D> void MWNode<D>::dealloc() {
    NOT_REACHED_ABORT;
}

/** Allocate the coefs vector. Only used by loose nodes. */
template <int D> void MWNode<D>::allocCoefs(int n_blocks, int block_size) {
    if (this->n_coefs != 0) MSG_ABORT("n_coefs should be zero");
    if (this->isAllocated()) MSG_ABORT("Coefs already allocated");
    if (not this->isLooseNode()) MSG_ABORT("Only loose nodes here!");

    this->n_coefs = n_blocks * block_size;
    this->coefs = new double[this->n_coefs];

    this->clearHasCoefs();
    this->setIsAllocated();
}

/** Deallocation of coefficients. Only used by loose nodes. */
template <int D> void MWNode<D>::freeCoefs() {
    if (not this->isLooseNode()) MSG_ABORT("Only loose nodes here!");

    if (this->coefs != nullptr) delete[] this->coefs;

    this->coefs = nullptr;
    this->n_coefs = 0;

    this->clearHasCoefs();
    this->clearIsAllocated();
}

template <int D> void MWNode<D>::printCoefs() const {
    if (not this->isAllocated()) MSG_ABORT("Node is not allocated");
    println(0, "\nMW coefs");
    int kp1_d = this->getKp1_d();
    for (int i = 0; i < this->n_coefs; i++) {
        println(0, this->coefs[i]);
        if (i % kp1_d == 0) println(0, "\n");
    }
}

template <int D> void MWNode<D>::getCoefs(Eigen::VectorXd &c) const {
    if (not this->isAllocated()) MSG_ABORT("Node is not allocated");
    if (not this->hasCoefs()) MSG_ABORT("Node has no coefs");
    if (this->n_coefs == 0) MSG_ABORT("ncoefs == 0");

    c = VectorXd::Map(this->coefs, this->n_coefs);
}

template <int D> void MWNode<D>::zeroCoefs() {
    if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");

    for (int i = 0; i < this->n_coefs; i++) { this->coefs[i] = 0.0; }
    this->zeroNorms();
    this->setHasCoefs();
}

template <int D> void MWNode<D>::setCoefBlock(int block, int block_size, const double *c) {
    if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
    for (int i = 0; i < block_size; i++) { this->coefs[block * block_size + i] = c[i]; }
}

template <int D> void MWNode<D>::addCoefBlock(int block, int block_size, const double *c) {
    if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
    for (int i = 0; i < block_size; i++) { this->coefs[block * block_size + i] += c[i]; }
}

template <int D> void MWNode<D>::zeroCoefBlock(int block, int block_size) {
    if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
    for (int i = 0; i < block_size; i++) { this->coefs[block * block_size + i] = 0.0; }
}

template <int D> void MWNode<D>::giveChildrenCoefs(bool overwrite) {
    assert(this->isBranchNode());
    if (not this->hasCoefs()) MSG_ABORT("No coefficients!");

    if (overwrite) {
        for (int i = 0; i < getTDim(); i++) getMWChild(i).zeroCoefs();
    }

    // coeff of child should be have been allocated already here
    int stride = getMWChild(0).getNCoefs();
    double *inp = getCoefs();
    double *out = getMWChild(0).getCoefs();
    bool readOnlyScaling = false;
    if (this->isGenNode()) readOnlyScaling = true;

    auto &tree = getMWTree();
    tree_utils::mw_transform(tree, inp, out, readOnlyScaling, stride, overwrite);

    for (int i = 0; i < getTDim(); i++) {
        getMWChild(i).setHasCoefs();
        getMWChild(i).calcNorms(); // should need to compute only scaling norms
    }
}

/** Takes the scaling coefficients of the children and stores them consecutively
 * in the  given vector. */
template <int D> void MWNode<D>::copyCoefsFromChildren() {
    int kp1_d = this->getKp1_d();
    int nChildren = this->getTDim();
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        MWNode<D> &child = getMWChild(cIdx);
        if (not child.hasCoefs()) MSG_ABORT("Child has no coefs");
        setCoefBlock(cIdx, kp1_d, child.getCoefs());
    }
}

/** Generates children nodes in a thread safe manner if (*this) is a leaf node */
template <int D> void MWNode<D>::threadSafeGenChildren() {
    // set_node_lock caused segmentation errors seldom and randomly -> avoided

    // we make many copies of the genChildren calls
    // the chance that 2 different nodes are treated in the same copy is small and therefore
    // the different critical sections will not interfer too much
    if (this->isLeafNode() or this->lockX != 0) {
        this->lockX = 1;

        int nCritical = 16;
        switch (this->serialIx % nCritical) {
            case 0:
#pragma omp critical(g0)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 1:
#pragma omp critical(g1)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 2:
#pragma omp critical(g2)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 3:
#pragma omp critical(g3)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 4:
#pragma omp critical(g4)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 5:
#pragma omp critical(g5)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 6:
#pragma omp critical(g6)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 7:
#pragma omp critical(g7)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 8:
#pragma omp critical(g8)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 9:
#pragma omp critical(g9)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 10:
#pragma omp critical(g10)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 11:
#pragma omp critical(g11)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 12:
#pragma omp critical(g12)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 13:
#pragma omp critical(g13)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 14:
#pragma omp critical(g14)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            case 15:
#pragma omp critical(g15)
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
                break;
            default:
#pragma omp critical
                if (isLeafNode()) {
                    genChildren();
                    giveChildrenCoefs();
                }
        }
        this->lockX = 0;
    }
}

/** Coefficient-Value transform
 *
 * This routine transforms the scaling coefficients of the node to the
 * function values in the corresponding quadrature roots (of its children).
 * Input parameter = forward: coef->value.
 * Input parameter = backward: value->coef.
 *
 * NOTE: this routine assumes a 0/1 (scaling on children 0 and 1)
 *       representation, in oppose to s/d (scaling and wavelet). */
template <int D> void MWNode<D>::cvTransform(int operation) {
    int kp1 = this->getKp1();
    int kp1_dm1 = math_utils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    int nCoefs = this->getTDim() * kp1_d;

    auto sb = this->getMWTree().getMRA().getScalingBasis();
    const MatrixXd &S = sb.getCVMap(operation);
    double o_vec[nCoefs];
    double *out_vec = o_vec;
    double *in_vec = this->coefs;

    for (int i = 0; i < D; i++) {
        for (int t = 0; t < this->getTDim(); t++) {
            double *out = out_vec + t * kp1_d;
            double *in = in_vec + t * kp1_d;
            math_utils::apply_filter(out, in, S, kp1, kp1_dm1, 0.0);
        }
        double *tmp = in_vec;
        in_vec = out_vec;
        out_vec = tmp;
    }

    const auto sf = this->getMWTree().getMRA().getWorldBox().getScalingFactor();
    double sf_prod = 1.0;
    for (const auto &s : sf) sf_prod *= s;
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
template<int D>
void MWNode<D>::cvTransform(int operation) {
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

/** Multiwavelet transform: fast version
  *
  * Application of the filters on one node to pass from a 0/1 (scaling
  * on children 0 and 1) representation to an s/d (scaling and
  * wavelet) representation. Bit manipulation is used in order to
  * determine the correct filters and whether to apply them or just
  * pass to the next couple of indexes. The starting coefficients are
  * preserved until the application is terminated, then they are
  * overwritten. With minor modifications this code can also be used
  * for the inverse mw transform (just use the transpose filters) or
  * for the application of an operator (using A, B, C and T parts of an
  * operator instead of G1, G0, H1, H0). This is the version where the
  * three directions are operated one after the other. Although this
  * is formally faster than the other algorithm, the separation of the
  * three dimensions prevent the possibility to use the norm of the
  * operator in order to discard a priori negligible contributions.

  * Luca Frediani, August 2006
  * C++ version: Jonas Juselius, September 2009 */
template <int D> void MWNode<D>::mwTransform(int operation) {
    int kp1 = this->getKp1();
    int kp1_dm1 = math_utils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    int nCoefs = this->getTDim() * kp1_d;
    const MWFilter &filter = getMWTree().getMRA().getFilter();
    double overwrite = 0.0;

    double o_vec[nCoefs];
    double *out_vec = o_vec;
    double *in_vec = this->coefs;

    for (int i = 0; i < D; i++) {
        int mask = 1 << i;
        for (int gt = 0; gt < this->getTDim(); gt++) {
            double *out = out_vec + gt * kp1_d;
            for (int ft = 0; ft < this->getTDim(); ft++) {
                /* Operate in direction i only if the bits along other
                 * directions are identical. The bit of the direction we
                 * operate on determines the appropriate filter/operator */
                if ((gt | mask) == (ft | mask)) {
                    double *in = in_vec + ft * kp1_d;
                    int fIdx = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                    const MatrixXd &oper = filter.getSubFilter(fIdx, operation);
                    math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                    overwrite = 1.0;
                }
            }
            overwrite = 0.0;
        }
        double *tmp = in_vec;
        in_vec = out_vec;
        out_vec = tmp;
    }
    if (IS_ODD(D)) {
        for (int i = 0; i < nCoefs; i++) { this->coefs[i] = in_vec[i]; }
    }
}

/** Set all norms to Undefined. */
template <int D> void MWNode<D>::clearNorms() {
    this->squareNorm = -1.0;
    for (int i = 0; i < this->getTDim(); i++) { this->componentNorms[i] = -1.0; }
}

/** Set all norms to zero. */
template <int D> void MWNode<D>::zeroNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) { this->componentNorms[i] = 0.0; }
}

/** Calculate and store square norm and component norms, if allocated. */
template <int D> void MWNode<D>::calcNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) {
        double norm_i = calcComponentNorm(i);
        this->componentNorms[i] = norm_i;
        this->squareNorm += norm_i * norm_i;
    }
}

/** Calculate and return the squared scaling norm. */
template <int D> double MWNode<D>::getScalingNorm() const {
    double sNorm = this->getComponentNorm(0);
    if (sNorm >= 0.0) {
        return sNorm * sNorm;
    } else {
        return -1.0;
    }
}

/** Calculate and return the squared wavelet norm. */
template <int D> double MWNode<D>::getWaveletNorm() const {
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

/** Calculate the norm of one component (NOT the squared norm!). */
template <int D> double MWNode<D>::calcComponentNorm(int i) const {
    assert(i == 0 or (not this->isGenNode())); // GenNodes have no wavelets coefficients
    assert(this->isAllocated());
    assert(this->hasCoefs());

    const double *c = this->getCoefs();
    int size = this->getKp1_d();
    int start = i * size;

    double sq_norm = 0.0;
#ifdef HAVE_BLAS
    sq_norm = cblas_ddot(size, &c[start], 1, &c[start], 1);
#else
    for (int i = start; i < start + size; i++) { sq_norm += c[i] * c[i]; }
#endif
    return std::sqrt(sq_norm);
}

/** Update the coefficients of the node by a mw transform of the scaling
 * coefficients of the children. Option to overwrite or add up existing
 * coefficients. */
template <int D> void MWNode<D>::reCompress() {
    if (this->isBranchNode()) {
        if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
        copyCoefsFromChildren();
        mwTransform(Compression);
        this->setHasCoefs();
        this->calcNorms();
    }
}

/** Recurse down until an EndNode is found, and then crop children with
 * too high precision. */
template <int D> bool MWNode<D>::crop(double prec, double splitFac, bool absPrec) {
    if (this->isEndNode()) {
        return true;
    } else {
        for (int i = 0; i < this->getTDim(); i++) {
            MWNode<D> &child = *this->children[i];
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

template <int D> void MWNode<D>::createChildren() {
    if (this->isBranchNode()) MSG_ABORT("Node already has children");
    this->getMWTree().getNodeAllocator().allocChildren(*this);
    this->setIsBranchNode();
}

template <int D> void MWNode<D>::createChildrenNoCoeff() {
    if (this->isBranchNode()) MSG_ABORT("Node already has children");
    this->getMWTree().getNodeAllocator().allocChildrenNoCoeff(*this);
    this->setIsBranchNode();
}

template <int D> void MWNode<D>::genChildren() {
    NOT_REACHED_ABORT;
}

/** Recursive deallocation of children and all their descendants.
 * Leaves node as LeafNode and children[] as null pointer. */
template <int D> void MWNode<D>::deleteChildren() {
    if (this->isLeafNode()) return;
    for (int cIdx = 0; cIdx < getTDim(); cIdx++) {
        if (this->children[cIdx] != nullptr) {
            MWNode<D> &child = getMWChild(cIdx);
            child.deleteChildren();
            child.dealloc();
        }
        this->children[cIdx] = nullptr;
    }
    this->childSerialIx = -1;
    this->setIsLeafNode();
}

template <int D> void MWNode<D>::deleteGenerated() {
    if (this->isBranchNode()) {
        if (this->isEndNode()) {
            this->deleteChildren();
        } else {
            for (int cIdx = 0; cIdx < getTDim(); cIdx++) { this->getMWChild(cIdx).deleteGenerated(); }
        }
    }
}

template <int D> void MWNode<D>::getCenter(double *r) const {
    NOT_IMPLEMENTED_ABORT;
    //    assert(r != 0);
    //    double sFac = std::pow(2.0, -getScale());
    //    for (int d = 0; d < D; d++) {
    //        double l = (double) getTranslation()[d];
    //        r[d] = sFac*(l + 0.5);
    //    }
}

template <int D> void MWNode<D>::getBounds(double *lb, double *ub) const {
    const auto sf = getMWTree().getMRA().getWorldBox().getScalingFactor();
    int n = getScale();
    double p = std::pow(2.0, -n);
    const NodeIndex<D> &l = getNodeIndex();
    for (int i = 0; i < D; i++) {
        lb[i] = sf[i] * (p * l[i]);
        ub[i] = sf[i] * (p * (l[i] + 1));
    }
}

/** Routine to find the path along the tree.
 *
 * Given the translation indices at the final scale, computes the child m
 * to be followed at the current scale in oder to get to the requested
 * node at the final scale. The result is the index of the child needed.
 * The index is obtained by bit manipulation of of the translation indices. */
template <int D> int MWNode<D>::getChildIndex(const NodeIndex<D> &nIdx) const {
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

/** Routine to find the path along the tree.
 *
 * Given a point in space, determines which child should be followed
 * to get to the corresponding terminal node. */
template <int D> int MWNode<D>::getChildIndex(const Coord<D> &r) const {
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

template <int D> void MWNode<D>::getPrimitiveQuadPts(MatrixXd &pts) const {
    int kp1 = this->getKp1();
    pts = MatrixXd::Zero(kp1, D);

    getQuadratureCache(qc);
    const VectorXd &roots = qc.getRoots(kp1);

    double sFac = std::pow(2.0, -getScale());
    const NodeIndex<D> &l = getNodeIndex();
    for (int d = 0; d < D; d++) pts.col(d) = sFac * (roots.array() + static_cast<double>(l[d]));
}

template <int D> void MWNode<D>::getPrimitiveChildPts(MatrixXd &pts) const {
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

template <int D> void MWNode<D>::getExpandedQuadPts(Eigen::MatrixXd &pts) const {
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

template <int D> void MWNode<D>::getExpandedChildPts(MatrixXd &pts) const {
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

/** Const version of node retriever that NEVER generates.
 *
 * Recursive routine to find and return the node with a given NodeIndex.
 * This routine returns the appropriate ProjectedNode, or a NULL pointer if
 * the node does not exist, or if it is a GenNode. Recursion starts at at this
 * node and ASSUMES the requested node is in fact decending from this node. */
template <int D> const MWNode<D> *MWNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) const {
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

/** Node retriever that NEVER generates.
 *
 * Recursive routine to find and return the node with a given NodeIndex.
 * This routine returns the appropriate ProjectedNode, or a NULL pointer if
 * the node does not exist, or if it is a GenNode. Recursion starts at at this
 * node and ASSUMES the requested node is in fact decending from this node. */
template <int D> MWNode<D> *MWNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) {
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

template <int D> const MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const Coord<D> &r, int depth) const {
    if (getDepth() == depth or this->isEndNode()) { return this; }
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNodeOrEndNode(r, depth);
}

/** Node retriever that return requested ProjectedNode or EndNode.
 *
 * Recursive routine to find and return the node with a given NodeIndex.
 * This routine returns the appropriate ProjectedNode, or the EndNode on the
 * path to the requested node, and will never create or return GenNodes.
 * Recursion starts at at this node and ASSUMES the requested node is in fact
 * decending from this node. */
template <int D> MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const Coord<D> &r, int depth) {
    if (getDepth() == depth or this->isEndNode()) { return this; }
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNodeOrEndNode(r, depth);
}

template <int D> const MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) const {
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

template <int D> MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) {
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

/** Node retriever that ALWAYS returns the requested node.
 *
 * Recursive routine to find and return the node with a given NodeIndex.
 * This routine always returns the appropriate node, and will generate nodes
 * that does not exist. Recursion starts at this node and ASSUMES the
 * requested node is in fact decending from this node. */
template <int D> MWNode<D> *MWNode<D>::retrieveNode(const Coord<D> &r, int depth) {
    if (depth < 0) MSG_ABORT("Invalid argument");

    if (getDepth() == depth) { return this; }
    assert(hasCoord(r));
    threadSafeGenChildren();
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNode(r, depth);
}

/** Node retriever that ALWAYS returns the requested node, possibly without coefs.
 *
 * Recursive routine to find and return the node with a given NodeIndex. This
 * routine always returns the appropriate node, and will generate nodes that
 * does not exist. Recursion starts at this node and ASSUMES the requested
 * node is in fact decending from this node. */
template <int D> MWNode<D> *MWNode<D>::retrieveNode(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    threadSafeGenChildren();
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != nullptr);
    return this->children[cIdx]->retrieveNode(idx);
}

/** Gives the norm (absolute value) of the node at the given NodeIndex.
 *
 * Recursive routine to find the node with a given NodeIndex. When an EndNode is
 * found, do not generate any new node, but rather give the value of the norm
 * assuming the function is uniformly distributed within the node. */
template <int D> double MWNode<D>::getNodeNorm(const NodeIndex<D> &idx) const {
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

/** Test if a given coordinate is within the boundaries of the node. */
template <int D> bool MWNode<D>::hasCoord(const Coord<D> &r) const {
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
template <int D> bool MWNode<D>::isCompatible(const MWNode<D> &node) {
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

/** Test if the node is decending from a given NodeIndex, that is, if they have
 * overlapping support. */
template <int D> bool MWNode<D>::isAncestor(const NodeIndex<D> &idx) const {
    int relScale = idx.getScale() - getScale();
    if (relScale < 0) return false;
    const NodeIndex<D> &l = getNodeIndex();
    for (int d = 0; d < D; d++) {
        int reqTransl = idx[d] >> relScale;
        if (l[d] != reqTransl) return false;
    }
    return true;
}

template <int D> bool MWNode<D>::isDecendant(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> std::ostream &MWNode<D>::print(std::ostream &o) const {
    std::string flags = "       ";
    o << getNodeIndex();
    if (isRootNode()) { flags[0] = 'R'; }
    if (isEndNode()) { flags[1] = 'E'; }
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
    if (isAllocated()) { flags[4] = 'A'; }
    if (hasCoefs()) { flags[5] = 'C'; }
    o << " " << flags;
    o << " sqNorm=" << this->squareNorm;
    if (hasCoefs()) {
        o << " Coefs={";
        o << getCoefs()[0] << ", " << getCoefs()[getNCoefs() - 1] << "}";
    }
    return o;
}

/** @brief recursively set maxSquaredNorm and maxWSquareNorm of parent and descendants
 *
 * @details normalization is such that a constant function gives constant value,
 * i.e. *not* same normalization as a squareNorm
 */
template <int D> void MWNode<D>::setMaxSquareNorm() {
    auto n = this->getScale();
    this->maxWSquareNorm = calcScaledWSquareNorm();
    this->maxSquareNorm = calcScaledSquareNorm();

    if (not this->isEndNode()) {
        for (int i = 0; i < this->getTDim(); i++) {
            MWNode<D> &child = *this->children[i];
            child.setMaxSquareNorm();
            this->maxSquareNorm = std::max(this->maxSquareNorm, child.maxSquareNorm);
            this->maxWSquareNorm = std::max(this->maxWSquareNorm, child.maxWSquareNorm);
        }
    }
}
/** @brief recursively reset maxSquaredNorm and maxWSquareNorm of parent and descendants to value -1 */
template <int D> void MWNode<D>::resetMaxSquareNorm() {
    auto n = this->getScale();
    this->maxSquareNorm = -1.0;
    this->maxWSquareNorm = -1.0;
    if (not this->isEndNode()) {
        for (int i = 0; i < this->getTDim(); i++) {
            MWNode<D> &child = *this->children[i];
            child.resetMaxSquareNorm();
        }
    }
}

template class MWNode<1>;
template class MWNode<2>;
template class MWNode<3>;

} // namespace mrcpp
