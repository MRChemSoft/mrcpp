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

#include "FunctionNode.h"
#include "FunctionTree.h"
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
template <int D> double FunctionNode<D>::evalf(Coord<D> r) {
    if (not this->hasCoefs()) MSG_ERROR("Evaluating node without coefs");

    // The 1.0 appearing in the if tests comes from the period is always 1.0
    // from the point of view of this function.
    if (this->getMWTree().getRootBox().isPeriodic()) {
        periodic::coord_manipulation<D>(r, this->getMWTree().getRootBox().getPeriodic());
    }

    this->threadSafeGenChildren();
    int cIdx = this->getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return getFuncChild(cIdx).evalScaling(r);
}

template <int D> double FunctionNode<D>::evalScaling(const Coord<D> &r) const {
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

    double result = 0.0;
    //#pragma omp parallel for shared(fact) reduction(+:result) num_threads(mrcpp_get_num_threads())
    for (int i = 0; i < this->getKp1_d(); i++) {
        double temp = this->coefs[i];
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
template <int D> double FunctionNode<D>::integrate() const {
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
template <int D> double FunctionNode<D>::integrateLegendre() const {
    double n = (D * this->getScale()) / 2.0;
    double two_n = std::pow(2.0, -n);
    return two_n * this->getCoefs()[0];
}

/** Function integration, Interpolating basis.
 *
 * Integrates the function represented on the node on the full support of the
 * node. A bit more involved than in the Legendre basis, as is requires some
 * coupling of quadrature weights. */
template <int D> double FunctionNode<D>::integrateInterpolating() const {
    int qOrder = this->getKp1();
    getQuadratureCache(qc);
    const VectorXd &weights = qc.getWeights(qOrder);

    double sqWeights[qOrder];
    for (int i = 0; i < qOrder; i++) sqWeights[i] = std::sqrt(weights[i]);

    int kp1_p[D];
    for (int i = 0; i < D; i++) kp1_p[i] = math_utils::ipow(qOrder, i);

    VectorXd coefs;
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
    double sum = coefs.segment(0, this->getKp1_d()).sum();

    return two_n * sum;
}

template <int D> void FunctionNode<D>::setValues(const VectorXd &vec) {
    this->zeroCoefs();
    this->setCoefBlock(0, vec.size(), vec.data());
    this->cvTransform(Backward);
    this->mwTransform(Compression);
    this->setHasCoefs();
    this->calcNorms();
}

template <int D> void FunctionNode<D>::getValues(VectorXd &vec) {
    if (this->isGenNode()) {
        MWNode<D> copy(*this);
        vec = Eigen::VectorXd::Zero(copy.getNCoefs());
        copy.mwTransform(Reconstruction);
        copy.cvTransform(Forward);
        for (int i = 0; i < this->n_coefs; i++) vec(i) = copy.getCoefs()[i];
    } else {
        vec = VectorXd::Zero(this->n_coefs);
        this->mwTransform(Reconstruction);
        this->cvTransform(Forward);
        for (int i = 0; i < this->n_coefs; i++) vec(i) = this->coefs[i];
        this->cvTransform(Backward);
        this->mwTransform(Compression);
    }
}

/** get coefficients corresponding to absolute value of function
 *
 * Leaves the original coefficients unchanged. */
template <int D> void FunctionNode<D>::getAbsCoefs(double *absCoefs) {
    double *coefsTmp = this->coefs;
    for (int i = 0; i < this->n_coefs; i++) absCoefs[i] = coefsTmp[i]; // copy
    this->coefs = absCoefs;                                            // swap coefs
    this->mwTransform(Reconstruction);
    this->cvTransform(Forward);
    for (int i = 0; i < this->n_coefs; i++) this->coefs[i] = std::abs(this->coefs[i]);
    this->cvTransform(Backward);
    this->mwTransform(Compression);
    this->coefs = coefsTmp; // restore original array (same address)
}

template <int D> void FunctionNode<D>::createChildren(bool coefs) {
    MWNode<D>::createChildren(coefs);
    this->clearIsEndNode();
}

template <int D> void FunctionNode<D>::genChildren() {
    if (this->isBranchNode()) MSG_ABORT("Node already has children");
    this->getFuncTree().getGenNodeAllocator().allocChildren(*this, true, true);
    this->setIsBranchNode();
}

template <int D> void FunctionNode<D>::deleteChildren() {
    MWNode<D>::deleteChildren();
    this->setIsEndNode();
}

template <int D> void FunctionNode<D>::dealloc() {
    int sIdx = this->serialIx;
    this->serialIx = -1;
    this->parentSerialIx = -1;
    this->childSerialIx = -1;
    auto &ftree = this->getFuncTree();
    if (this->isGenNode()) {
        ftree.getGenNodeAllocator().deallocNodes(sIdx);
    } else {
        ftree.decrementNodeCount(this->getScale());
        ftree.getNodeAllocator().deallocNodes(sIdx);
    }
}

/** Update the coefficients of the node by a mw transform of the scaling
 * coefficients of the children. Option to overwrite or add up existing
 * coefficients. Specialized for D=3 below. */
template <int D> void FunctionNode<D>::reCompress() {
    MWNode<D>::reCompress();
}

template <> void FunctionNode<3>::reCompress() {
    if (this->isBranchNode()) {
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
 * coefficient vectors. Assumes the nodes have identical support. */
template <int D> double dot_scaling(const FunctionNode<D> &bra, const FunctionNode<D> &ket) {
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

/** Inner product of the functions represented by the wavelet basis of the nodes.
 *
 * Integrates the product of the functions represented by the wavelet basis on
 * the node on the full support of the nodes. The wavelet basis is fully
 * orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support. */
template <int D> double dot_wavelet(const FunctionNode<D> &bra, const FunctionNode<D> &ket) {
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

template double dot_scaling(const FunctionNode<1> &bra, const FunctionNode<1> &ket);
template double dot_scaling(const FunctionNode<2> &bra, const FunctionNode<2> &ket);
template double dot_scaling(const FunctionNode<3> &bra, const FunctionNode<3> &ket);
template double dot_wavelet(const FunctionNode<1> &bra, const FunctionNode<1> &ket);
template double dot_wavelet(const FunctionNode<2> &bra, const FunctionNode<2> &ket);
template double dot_wavelet(const FunctionNode<3> &bra, const FunctionNode<3> &ket);

template class FunctionNode<1>;
template class FunctionNode<2>;
template class FunctionNode<3>;

} // namespace mrcpp
