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

#include "FunctionNode.h"
#include "FunctionTree.h"
#include "NodeAllocator.h"
#include "core/QuadratureCache.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"
#include "utils/tree_utils.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace Eigen;

namespace mrcpp {

/** Function evaluation.
 * Evaluate all polynomials defined on the node. */
template <int D, typename T> T FunctionNode<D, T>::evalf(Coord<D> r) {
    if (not this->hasCoefs()) MSG_ERROR("Evaluating node without coefs");

    // The 1.0 appearing in the if tests comes from the period is always 1.0
    // from the point of view of this function.
    if (this->getMWTree().getRootBox().isPeriodic()) { periodic::coord_manipulation<D>(r, this->getMWTree().getRootBox().getPeriodic()); }

    this->threadSafeGenChildren();
    int cIdx = this->getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return getFuncChild(cIdx).evalScaling(r);
}

template <int D, typename T> T FunctionNode<D, T>::evalScaling(const Coord<D> &r) const {
    if (not this->hasCoefs()) MSG_ERROR("Evaluating node without coefs");

    double arg[D];
    double n_factor = std::pow(2.0, this->getScale());
    const NodeIndex<D> &l_factor = this->getNodeIndex();
    for (int i = 0; i < D; i++) arg[i] = r[i] * n_factor - static_cast<double>(l_factor[i]);

    int fact[D + 1];
    for (int i = 0; i < D + 1; i++) fact[i] = math_utils::ipow(this->getKp1(), i);

    MatrixXd val(this->getKp1(), D);
    const ScalingBasis &basis = this->getMWTree().getMRA().getScalingBasis();
    basis.evalf(arg, val);

    T result = 0.0;
    //#pragma omp parallel for shared(fact) reduction(+:result) num_threads(mrcpp_get_num_threads())
    for (int i = 0; i < this->getKp1_d(); i++) {
        T temp = this->coefs[i];
        for (int j = 0; j < D; j++) {
            int k = (i % fact[j + 1]) / fact[j];
            temp *= val(k, j);
        }
        result += temp;
    }
    double n = (D * this->getScale()) / 2.0;
    double two_n = std::pow(2.0, n);
    return two_n * result;
}

/** Function integration.
 *
 * Wrapper for function integration, that requires different methods depending
 * on scaling type. Integrates the function represented on the node on the
 * full support of the node. */
template <int D, typename T> T FunctionNode<D, T>::integrate() const {
    if (not this->hasCoefs()) { return 0.0; }
    switch (this->getScalingType()) {
        case Legendre:
            return integrateLegendre();
            break;
        case Interpol:
            return integrateInterpolating();
            break;
        default:
            MSG_ABORT("Invalid scalingType");
    }
}

/** Function integration, Legendre basis.
 *
 * Integrates the function represented on the node on the full support of the
 * node. The Legendre basis is particularly easy to integrate, as the work is
 * already done when calculating its coefficients. The coefficients of the
 * node is defined as the projection integral
 *          s_i = int f(x)phi_i(x)dx
 * and since the first Legendre function is the constant 1, the first
 * coefficient is simply the integral of f(x). */
template <int D, typename T> T FunctionNode<D, T>::integrateLegendre() const {
    double n = (D * this->getScale()) / 2.0;
    double two_n = std::pow(2.0, -n);
    return two_n * this->getCoefs()[0];
}

/** Function integration, Interpolating basis.
 *
 * Integrates the function represented on the node on the full support of the
 * node. A bit more involved than in the Legendre basis, as is requires some
 * coupling of quadrature weights. */
template <int D, typename T> T FunctionNode<D, T>::integrateInterpolating() const {
    int qOrder = this->getKp1();
    getQuadratureCache(qc);
    const VectorXd &weights = qc.getWeights(qOrder);
    std::vector<double> sqWeights(qOrder);
    for (int i = 0; i < qOrder; i++) sqWeights[i] = std::sqrt(weights[i]);

    int kp1_p[D];
    for (int i = 0; i < D; i++) kp1_p[i] = math_utils::ipow(qOrder, i);

    Eigen::Matrix<T, Eigen::Dynamic, 1> coefs;
    this->getCoefs(coefs);
    for (int p = 0; p < D; p++) {

        int n = 0;
        for (int i = 0; i < kp1_p[D - p - 1]; i++) {
            for (int j = 0; j < qOrder; j++) {
                for (int k = 0; k < kp1_p[p]; k++) {
                    coefs(n) *= sqWeights[j];
                    n++;
                }
            }
        }
    }
    double n = (D * this->getScale()) / 2.0;
    double two_n = std::pow(2.0, -n);
    T sum = coefs.segment(0, this->getKp1_d()).sum();

    return two_n * sum;
}

/** Function integration, Interpolating basis.
 *
 * Integrates the function represented on the node on the full support of the
 * node. A bit more involved than in the Legendre basis, as is requires some
 * coupling of quadrature weights. */
template <int D, typename T> T FunctionNode<D, T>::integrateValues() const {
    int qOrder = this->getKp1();
    getQuadratureCache(qc);
    const VectorXd &weights = qc.getWeights(qOrder);
    Eigen::Matrix<T, Eigen::Dynamic, 1> coefs;
    this->getCoefs(coefs);
    int ncoefs = coefs.size();
    int ncoefChild = ncoefs / (1 << D);
    std::vector<T> cc(ncoefChild);
    // factorize out the children
    for (int i = 0; i < ncoefChild; i++) cc[i] = coefs[i];
    for (int j = 1; j < (1 << D); j++)
        for (int i = 0; i < ncoefChild; i++) cc[i] += coefs[j * ncoefChild + i];

    int nc = 0;
    T sum = 0.0;
    if (D > 3)
        MSG_ABORT("Not Implemented")
    else if (D == 3) {
        for (int i = 0; i < qOrder; i++) {
            T sumj = 0.0;
            for (int j = 0; j < qOrder; j++) {
                T sumk = 0.0;
                for (int k = 0; k < qOrder; k++) sumk += cc[nc++] * weights[k];
                sumj += sumk * weights[j];
            }
            sum += sumj * weights[i];
        }
    } else if (D == 2) {
        for (int j = 0; j < qOrder; j++) {
            T sumk = 0.0;
            for (int k = 0; k < qOrder; k++) sumk += cc[nc++] * weights[k];
            sum += sumk * weights[j];
        }
    } else if (D == 1)
        for (int k = 0; k < qOrder; k++) sum += cc[nc++] * weights[k];

    int n = D * (this->getScale() + 1); // NB: one extra scale
    int two_n = (1 << abs(n));          // 2**n;
    if (n > 0)
        sum /= two_n;
    else
        sum *= two_n;
    return sum;
}

template <int D, typename T> void FunctionNode<D, T>::setValues(const Matrix<T, Eigen::Dynamic, 1> &vec) {
    this->zeroCoefs();
    this->setCoefBlock(0, vec.size(), vec.data());
    this->cvTransform(Backward);
    this->mwTransform(Compression);
    this->setHasCoefs();
    this->calcNorms();
}

template <int D, typename T> void FunctionNode<D, T>::getValues(Matrix<T, Eigen::Dynamic, 1> &vec) {
    if (this->isGenNode()) {
        MWNode<D, T> copy(*this);
        vec = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(copy.getNCoefs());
        copy.mwTransform(Reconstruction);
        copy.cvTransform(Forward);
        for (int i = 0; i < this->n_coefs; i++) vec(i) = copy.getCoefs()[i];
    } else {
        vec = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(this->n_coefs);
        this->mwTransform(Reconstruction);
        this->cvTransform(Forward);
        for (int i = 0; i < this->n_coefs; i++) vec(i) = this->coefs[i];
        this->cvTransform(Backward);
        this->mwTransform(Compression);
    }
}

/** get coefficients corresponding to absolute value of function
 *
 * Leaves the original coefficients unchanged.
 * Note that we mus use T and not double, even if the norms are double, because
 * the transforms expect T types.
 */
template <int D, typename T> void FunctionNode<D, T>::getAbsCoefs(T *absCoefs) {
    T *coefsTmp = this->coefs;
    for (int i = 0; i < this->n_coefs; i++) absCoefs[i] = coefsTmp[i]; // copy
    this->coefs = absCoefs;                                            // swap coefs
    this->mwTransform(Reconstruction);
    this->cvTransform(Forward);
    for (int i = 0; i < this->n_coefs; i++) this->coefs[i] = std::norm(this->coefs[i]);
    this->cvTransform(Backward);
    this->mwTransform(Compression);
    this->coefs = coefsTmp; // restore original array (same address)
}

template <int D, typename T> void FunctionNode<D, T>::createChildren(bool coefs) {
    if (this->isBranchNode()) MSG_ABORT("Node already has children");
    auto &allocator = this->getFuncTree().getNodeAllocator();

    int nChildren = this->getTDim();
    int sIdx = allocator.alloc(nChildren, coefs);

    auto n_coefs = allocator.getNCoefs();
    auto *coefs_p = (coefs) ? allocator.getCoef_p(sIdx) : nullptr;
    auto *child_p = allocator.getNode_p(sIdx);

    this->childSerialIx = sIdx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        // construct into allocator memory
        new (child_p) FunctionNode<D, T>(this, cIdx);
        this->children[cIdx] = child_p;

        child_p->serialIx = sIdx;
        child_p->parentSerialIx = this->serialIx;
        child_p->childSerialIx = -1;

        child_p->n_coefs = n_coefs;
        child_p->coefs = coefs_p;
        if (coefs) child_p->setIsAllocated();

        child_p->setIsLeafNode();
        child_p->setIsEndNode();
        child_p->clearHasCoefs();

        this->getMWTree().incrementNodeCount(child_p->getScale());
        sIdx++;
        child_p++;
        if (coefs) coefs_p += n_coefs;
    }
    this->setIsBranchNode();
    this->clearIsEndNode();
}

template <int D, typename T> void FunctionNode<D, T>::genChildren() {
    if (this->isBranchNode()) MSG_ABORT("Node already has children");
    auto &allocator = this->getFuncTree().getGenNodeAllocator();

    int nChildren = this->getTDim();
    int sIdx = allocator.alloc(nChildren);

    auto n_coefs = allocator.getNCoefs();
    auto *coefs_p = allocator.getCoef_p(sIdx);
    auto *child_p = allocator.getNode_p(sIdx);

    this->childSerialIx = sIdx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        // construct into allocator memory
        new (child_p) FunctionNode<D, T>(this, cIdx);
        this->children[cIdx] = child_p;

        child_p->serialIx = sIdx;
        child_p->parentSerialIx = (this->isGenNode()) ? this->serialIx : -1;
        child_p->childSerialIx = -1;

        child_p->n_coefs = n_coefs;
        child_p->coefs = coefs_p;
        child_p->setIsAllocated();

        child_p->setIsLeafNode();
        child_p->setIsGenNode();
        child_p->clearHasCoefs();
        child_p->clearIsEndNode();

        sIdx++;
        child_p++;
        coefs_p += n_coefs;
    }
    this->setIsBranchNode();
}

template <int D, typename T> void FunctionNode<D, T>::genParent() {
    if (this->parent != nullptr) MSG_ABORT("Node is not an orphan");

    auto &allocator = this->getFuncTree().getNodeAllocator();
    int sIdx = allocator.alloc(1, true);

    auto n_coefs = allocator.getNCoefs();
    auto *coefs_p = allocator.getCoef_p(sIdx);
    auto *parent_p = allocator.getNode_p(sIdx);

    this->parentSerialIx = sIdx;

    // construct into allocator memory
    new (parent_p) FunctionNode<D, T>(this->tree, this->getNodeIndex().parent());

    this->parent = parent_p;

    for (int cIdx = 0; cIdx < this->getTDim(); cIdx++) parent_p->children[cIdx] = this;
    parent_p->serialIx = sIdx;
    parent_p->parentSerialIx = -1;
    parent_p->childSerialIx = this->serialIx;

    parent_p->n_coefs = n_coefs;
    parent_p->coefs = coefs_p;

    parent_p->setIsBranchNode();
    parent_p->setIsAllocated();
    parent_p->clearHasCoefs();

    this->getMWTree().incrementNodeCount(parent_p->getScale());
}

template <int D, typename T> void FunctionNode<D, T>::deleteChildren() {
    MWNode<D, T>::deleteChildren();
    this->setIsEndNode();
}

template <int D, typename T> void FunctionNode<D, T>::dealloc() {
    int sIdx = this->serialIx;
    this->serialIx = -1;
    this->parentSerialIx = -1;
    this->childSerialIx = -1;
    auto &ftree = this->getFuncTree();
    if (this->isGenNode()) {
        ftree.getGenNodeAllocator().dealloc(sIdx);
        // for GenNodes we clean unused chunks carefully, as they can become
        // very large and occupy space long after used.
        ftree.getGenNodeAllocator().deleteUnusedChunks();
    } else {
        ftree.decrementNodeCount(this->getScale());
        ftree.getNodeAllocator().dealloc(sIdx);
    }
}

/** Update the coefficients of the node by a mw transform of the scaling
 * coefficients of the children. Option to overwrite or add up existing
 * coefficients. Specialized for D=3 below. */
template <int D, typename T> void FunctionNode<D, T>::reCompress() {
    MWNode<D, T>::reCompress();
}

template <> void FunctionNode<3>::reCompress() {
    if (this->getDepth() < 0) {
        // This happens for negative scale pbc operators
        MWNode<3>::reCompress();
    } else if (this->isBranchNode()) {
        if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
        // can write directly from children coeff into parent coeff
        int stride = this->getMWChild(0).getNCoefs();
        double *inp = this->getMWChild(0).getCoefs();
        double *out = this->coefs;

        assert(inp + 7 * stride == this->getMWChild(7).getCoefs());

        auto &tree = getMWTree();
        tree_utils::mw_transform_back(tree, inp, out, stride);
        this->setHasCoefs();
        this->calcNorms();
    }
}

/** Inner product of the functions represented by the scaling basis of the nodes.
 *
 * Integrates the product of the functions represented by the scaling basis on
 * the node on the full support of the nodes. The scaling basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * NB: will take conjugate of bra in case of complex values.
 */
template <int D> double dot_scaling(const FunctionNode<D, double> &bra, const FunctionNode<D, double> &ket) {
    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const double *a = bra.getCoefs();
    const double *b = ket.getCoefs();

    int size = bra.getKp1_d();
#ifdef HAVE_BLAS
    return cblas_ddot(size, a, 1, b, 1);
#else
    double result = 0.0;
    for (int i = 0; i < size; i++) result += a[i] * b[i];
    return result;
#endif
}

/** Inner product of the functions represented by the scaling basis of the nodes.
 *
 * Integrates the product of the functions represented by the scaling basis on
 * the node on the full support of the nodes. The scaling basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * NB: will take conjugate of bra in case of complex values.
 */
template <int D> ComplexDouble dot_scaling(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, ComplexDouble> &ket) {
    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const ComplexDouble *a = bra.getCoefs();
    const ComplexDouble *b = ket.getCoefs();

    int size = bra.getKp1_d();
    ComplexDouble result = 0.0;
    // note that bra is conjugated by default
    if (bra.getMWTree().conjugate()) {
        if (ket.getMWTree().conjugate()) {
            for (int i = 0; i < size; i++) result += a[i] * std::conj(b[i]);
        } else {
            for (int i = 0; i < size; i++) result += a[i] * b[i];
        }
    } else {
        if (ket.getMWTree().conjugate()) {
            for (int i = 0; i < size; i++) result += std::conj(a[i]) * std::conj(b[i]);
        } else {
            for (int i = 0; i < size; i++) result += std::conj(a[i]) * b[i];
        }
    }
    return result;
}

/** Inner product of the functions represented by the scaling basis of the nodes.
 *
 * Integrates the product of the functions represented by the scaling basis on
 * the node on the full support of the nodes. The scaling basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * NB: will take conjugate of bra in case of complex values.
 */
template <int D> ComplexDouble dot_scaling(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, double> &ket) {
    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const ComplexDouble *a = bra.getCoefs();
    const double *b = ket.getCoefs();

    int size = bra.getKp1_d();
    ComplexDouble result = 0.0;
    // note that bra is conjugated by default
    if (bra.getMWTree().conjugate()) {
        for (int i = 0; i < size; i++) result += a[i] * b[i];
    } else {
        for (int i = 0; i < size; i++) result += std::conj(a[i]) * b[i];
    }
    return result;
}

/** Inner product of the functions represented by the scaling basis of the nodes.
 *
 * Integrates the product of the functions represented by the scaling basis on
 * the node on the full support of the nodes. The scaling basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * NB: will take conjugate of bra in case of complex values.
 */
template <int D> ComplexDouble dot_scaling(const FunctionNode<D, double> &bra, const FunctionNode<D, ComplexDouble> &ket) {
    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const double *a = bra.getCoefs();
    const ComplexDouble *b = ket.getCoefs();

    int size = bra.getKp1_d();
    ComplexDouble result = 0.0;
    // note that bra is conjugated by default
    if (ket.getMWTree().conjugate()) {
        for (int i = 0; i < size; i++) result += a[i] * std::conj(b[i]);
    } else {
        for (int i = 0; i < size; i++) result += a[i] * b[i];
    }
    return result;
}

/** Inner product of the functions represented by the wavelet basis of the nodes.
 *
 * Integrates the product of the functions represented by the wavelet basis on
 * the node on the full support of the nodes. The wavelet basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * NB: will take conjugate of bra in case of complex values.
 */
template <int D> double dot_wavelet(const FunctionNode<D, double> &bra, const FunctionNode<D, double> &ket) {
    if (bra.isGenNode() or ket.isGenNode()) return 0.0;

    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const double *a = bra.getCoefs();
    const double *b = ket.getCoefs();

    int start = bra.getKp1_d();
    int size = (bra.getTDim() - 1) * start;
#ifdef HAVE_BLAS
    return cblas_ddot(size, &a[start], 1, &b[start], 1);
#else
    double result = 0.0;
    for (int i = 0; i < size; i++) result += a[start + i] * b[start + i];
    return result;
#endif
}

/** Inner product of the functions represented by the wavelet basis of the nodes.
 *
 * Integrates the product of the functions represented by the wavelet basis on
 * the node on the full support of the nodes. The wavelet basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * NB: will take conjugate of bra in case of complex values.
 */
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, ComplexDouble> &ket) {
    if (bra.isGenNode() or ket.isGenNode()) return 0.0;

    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const ComplexDouble *a = bra.getCoefs();
    const ComplexDouble *b = ket.getCoefs();

    int start = bra.getKp1_d();
    int size = (bra.getTDim() - 1) * start;
    ComplexDouble result = 0.0;
    if (bra.getMWTree().conjugate()) {
        if (ket.getMWTree().conjugate()) {
            for (int i = 0; i < size; i++) result += a[start + i] * std::conj(b[start + i]);
        } else {
            for (int i = 0; i < size; i++) result += a[start + i] * b[start + i];
        }
    } else {
        if (ket.getMWTree().conjugate()) {
            for (int i = 0; i < size; i++) result += std::conj(a[start + i]) * std::conj(b[start + i]);
        } else {
            for (int i = 0; i < size; i++) result += std::conj(a[start + i]) * b[start + i];
        }
    }
    return result;
}

/** Inner product of the functions represented by the wavelet basis of the nodes.
 *
 * Integrates the product of the functions represented by the wavelet basis on
 * the node on the full support of the nodes. The wavelet basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * NB: will take conjugate of bra in case of complex values.
 */
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, double> &ket) {
    if (bra.isGenNode() or ket.isGenNode()) return 0.0;

    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const ComplexDouble *a = bra.getCoefs();
    const double *b = ket.getCoefs();

    int start = bra.getKp1_d();
    int size = (bra.getTDim() - 1) * start;
    ComplexDouble result = 0.0;
    if (bra.getMWTree().conjugate()) {
        for (int i = 0; i < size; i++) result += a[start + i] * b[start + i];
    } else {
        for (int i = 0; i < size; i++) result += std::conj(a[start + i]) * b[start + i];
    }
    return result;
}

/** Inner product of the functions represented by the wavelet basis of the nodes.
 *
 * Integrates the product of the functions represented by the wavelet basis on
 * the node on the full support of the nodes. The wavelet basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * NB: will take conjugate of bra in case of complex values.
 */
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, double> &bra, const FunctionNode<D, ComplexDouble> &ket) {
    if (bra.isGenNode() or ket.isGenNode()) return 0.0;

    assert(bra.hasCoefs());
    assert(ket.hasCoefs());

    const double *a = bra.getCoefs();
    const ComplexDouble *b = ket.getCoefs();

    int start = bra.getKp1_d();
    int size = (bra.getTDim() - 1) * start;
    ComplexDouble result = 0.0;
    if (ket.getMWTree().conjugate()) {
        for (int i = 0; i < size; i++) result += a[start + i] * std::conj(b[start + i]);
    } else {
        for (int i = 0; i < size; i++) result += a[start + i] * b[start + i];
    }
    return result;
}

template double dot_scaling(const FunctionNode<1, double> &bra, const FunctionNode<1, double> &ket);
template double dot_scaling(const FunctionNode<2, double> &bra, const FunctionNode<2, double> &ket);
template double dot_scaling(const FunctionNode<3, double> &bra, const FunctionNode<3, double> &ket);
template double dot_wavelet(const FunctionNode<1, double> &bra, const FunctionNode<1, double> &ket);
template double dot_wavelet(const FunctionNode<2, double> &bra, const FunctionNode<2, double> &ket);
template double dot_wavelet(const FunctionNode<3, double> &bra, const FunctionNode<3, double> &ket);

template class FunctionNode<1, double>;
template class FunctionNode<2, double>;
template class FunctionNode<3, double>;

template class FunctionNode<1, ComplexDouble>;
template class FunctionNode<2, ComplexDouble>;
template class FunctionNode<3, ComplexDouble>;

template ComplexDouble dot_scaling(const FunctionNode<1, ComplexDouble> &bra, const FunctionNode<1, ComplexDouble> &ket);
template ComplexDouble dot_scaling(const FunctionNode<2, ComplexDouble> &bra, const FunctionNode<2, ComplexDouble> &ket);
template ComplexDouble dot_scaling(const FunctionNode<3, ComplexDouble> &bra, const FunctionNode<3, ComplexDouble> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<1, ComplexDouble> &bra, const FunctionNode<1, ComplexDouble> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<2, ComplexDouble> &bra, const FunctionNode<2, ComplexDouble> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<3, ComplexDouble> &bra, const FunctionNode<3, ComplexDouble> &ket);

template ComplexDouble dot_scaling(const FunctionNode<1, double> &bra, const FunctionNode<1, ComplexDouble> &ket);
template ComplexDouble dot_scaling(const FunctionNode<2, double> &bra, const FunctionNode<2, ComplexDouble> &ket);
template ComplexDouble dot_scaling(const FunctionNode<3, double> &bra, const FunctionNode<3, ComplexDouble> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<1, double> &bra, const FunctionNode<1, ComplexDouble> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<2, double> &bra, const FunctionNode<2, ComplexDouble> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<3, double> &bra, const FunctionNode<3, ComplexDouble> &ket);

template ComplexDouble dot_scaling(const FunctionNode<1, ComplexDouble> &bra, const FunctionNode<1, double> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<1, ComplexDouble> &bra, const FunctionNode<1, double> &ket);
template ComplexDouble dot_scaling(const FunctionNode<2, ComplexDouble> &bra, const FunctionNode<2, double> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<2, ComplexDouble> &bra, const FunctionNode<2, double> &ket);
template ComplexDouble dot_scaling(const FunctionNode<3, ComplexDouble> &bra, const FunctionNode<3, double> &ket);
template ComplexDouble dot_wavelet(const FunctionNode<3, ComplexDouble> &bra, const FunctionNode<3, double> &ket);

} // namespace mrcpp
