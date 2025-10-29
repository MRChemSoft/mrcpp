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
 * @file DerivativeCalculator.cpp
 * @brief Node-wise application of multiresolution **derivative operators**.
 *
 * @details
 * This module implements the computational kernels used to apply a
 * #mrcpp::DerivativeOperator to multiresolution coefficient trees.
 * The calculator works node-by-node and supports both:
 *
 * - **Zero-bandwidth (local) operators** — e.g. ABGV-00 type operators that
 *   act diagonally (per cell) in non-applied directions; handled by
 *   DerivativeCalculator::applyOperator_bw0().
 * - **Finite-bandwidth operators** — e.g. ABGV-55/PH/BS operators that couple
 *   nearest neighbors along the application direction; handled by
 *   DerivativeCalculator::applyOperator().
 *
 * The apply pipeline for each output (g) node:
 *  1. Build the **operator band** of input (f) nodes affected by the operator
 *     at the current depth (makeOperBand()).
 *  2. For each combination of tensor components \f$(f_t,g_t)\f$, gather the
 *     1D operator blocks from the pre-built #mrcpp::OperatorTree components.
 *  3. Perform the separated **tensor contraction**
 *     (tensorApplyOperComp()) across dimensions, using identity where the
 *     operator is not applied.
 *  4. Apply a **scaling normalization** based on the world-box scaling factor
 *     and the derivative order.
 *  5. Compute node norms for downstream thresholding and diagnostics.
 *
 * The class collects **per-thread timing** and **operator-usage statistics**
 * to aid profiling and load balancing.
 *
 * @note Scaling normalization:
 *       The derivative w.r.t. direction `applyDir` is normalized by the
 *       world-box scaling factor of that direction raised to the operator
 *       order. See the notes near the end of calcNode() overloads.
 */

#include "DerivativeCalculator.h"
#include "operators/DerivativeOperator.h"
#include "operators/OperatorState.h"
#include "trees/BandWidth.h"
#include "trees/FunctionTree.h"
#include "trees/OperatorNode.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using Eigen::MatrixXd;

namespace mrcpp {

/**
 * @brief Construct a derivative calculator.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type (e.g., double or ComplexDouble).
 * @param dir  Direction along which the derivative is applied (0-based, < D).
 * @param o    Derivative operator to apply.
 * @param f    Input function tree (source of coefficients).
 *
 * @throws Aborts if `dir` is outside \f$[0,D)\f$.
 *
 * @details
 * The constructor stores references to the operator and the input tree,
 * validates the application direction, and initializes per-thread timers.
 */
template <int D, typename T>
DerivativeCalculator<D, T>::DerivativeCalculator(int dir, DerivativeOperator<D> &o, FunctionTree<D, T> &f)
        : applyDir(dir)
        , fTree(&f)
        , oper(&o) {
    if (dir < 0 or dir >= D) MSG_ABORT("Invalid apply dir");
    initTimers();
}

/**
 * @brief Flush usage counters and print aggregate statistics on destruction.
 */
template <int D, typename T> DerivativeCalculator<D, T>::~DerivativeCalculator() {
    this->operStat.flushNodeCounters();
    println(10, this->operStat);
}

/**
 * @brief Initialize per-thread timers (band construction / calc / norms).
 */
template <int D, typename T> void DerivativeCalculator<D, T>::initTimers() {
    int nThreads = mrcpp_get_max_threads();
    for (int i = 0; i < nThreads; i++) {
        this->band_t.push_back(Timer(false));
        this->calc_t.push_back(Timer(false));
        this->norm_t.push_back(Timer(false));
    }
}

/**
 * @brief Clear local timer storage.
 */
template <int D, typename T> void DerivativeCalculator<D, T>::clearTimers() {
    this->band_t.clear();
    this->calc_t.clear();
    this->norm_t.clear();
}

/**
 * @brief Print per-thread timing statistics gathered during application.
 */
template <int D, typename T> void DerivativeCalculator<D, T>::printTimers() const {
    int oldprec = Printer::setPrecision(1);
    int nThreads = mrcpp_get_max_threads();
    printout(20, "\n\nthread ");
    for (int i = 0; i < nThreads; i++) printout(20, std::setw(9) << i);
    printout(20, "\nband     ");
    for (int i = 0; i < nThreads; i++) printout(20, this->band_t[i].elapsed() << "  ");
    printout(20, "\ncalc     ");
    for (int i = 0; i < nThreads; i++) printout(20, this->calc_t[i].elapsed() << "  ");
    printout(20, "\nnorm     ");
    for (int i = 0; i < nThreads; i++) printout(20, this->norm_t[i].elapsed() << "  ");
    printout(20, "\n\n");
    Printer::setPrecision(oldprec);
}

/**
 * @brief Apply a **local (zero-bandwidth)** derivative operator to a single node.
 *
 * @param[in]  inpNode Source node (input function).
 * @param[out] outNode Destination node (output after derivative).
 *
 * @details
 * Uses applyOperator_bw0() which assumes the operator couples only the same
 * spatial cell in non-applied directions. Identity is implicitly used in
 * directions other than `applyDir`.
 *
 * After the tensor contraction, multiplies by \f$(\text{scale-factor})^{-p}\f$
 * where `p = oper->getOrder()` to account for physical coordinate scaling,
 * then updates norms.
 */
template <int D, typename T> void DerivativeCalculator<D, T>::calcNode(MWNode<D, T> &inpNode, MWNode<D, T> &outNode) {
    // if (this->oper->getMaxBandWidth() > 1) MSG_ABORT("Only implemented for zero bw");
    outNode.zeroCoefs();
    int nComp = (1 << D);
    T tmpCoefs[outNode.getNCoefs()];
    OperatorState<D, T> os(outNode, tmpCoefs);

    os.setFNode(inpNode);
    os.setFIndex(inpNode.nodeIndex);
    for (int ft = 0; ft < nComp; ft++) {
        double fNorm = inpNode.getComponentNorm(ft);
        if (fNorm < MachineZero) { continue; } // early skip of empty components
        os.setFComponent(ft);
        for (int gt = 0; gt < nComp; gt++) {
            os.setGComponent(gt);
            applyOperator_bw0(os);
        }
    }
    // Coordinate scaling normalization: divide by (scaleFactor^order).
    const double scaling_factor = 1.0 / std::pow(outNode.getMWTree().getMRA().getWorldBox().getScalingFactor(this->applyDir), oper->getOrder());
    if (abs(scaling_factor - 1.0) > MachineZero) {
        for (int i = 0; i < outNode.getNCoefs(); i++) outNode.getCoefs()[i] *= scaling_factor;
    }
    outNode.calcNorms(); // norms are used by downstream screening
}

/**
 * @brief Apply a **finite-bandwidth** derivative operator to a single node.
 *
 * @param gNode Destination/output node (will be overwritten).
 *
 * @details
 *  1. Build the operator band (list of input nodes influencing `gNode`) along
 *     the apply direction (makeOperBand()).
 *  2. For each band node and tensor-component pair, gather operator slices
 *     from the #mrcpp::OperatorTree and perform the tensor product.
 *  3. Apply coordinate scaling normalization by dividing by
 *     \f$(\text{scaleFactor})^{\text{order}}\f$ in the applied direction.
 *  4. Update norms and timing statistics.
 */
template <int D, typename T> void DerivativeCalculator<D, T>::calcNode(MWNode<D, T> &gNode) {
    gNode.zeroCoefs();

    int nComp = (1 << D);
    T tmpCoefs[gNode.getNCoefs()];
    OperatorState<D, T> os(gNode, tmpCoefs);
    this->operStat.incrementGNodeCounters(gNode);

    // Build band of input nodes that affect gNode
    this->band_t[mrcpp_get_thread_num()].resume();
    std::vector<NodeIndex<D>> idx_band;
    MWNodeVector<D, T> fBand = makeOperBand(gNode, idx_band);
    this->band_t[mrcpp_get_thread_num()].stop();

    this->calc_t[mrcpp_get_thread_num()].resume();
    for (int n = 0; n < fBand.size(); n++) {
        MWNode<D, T> &fNode = *fBand[n];
        NodeIndex<D> &fIdx = idx_band[n];
        os.setFNode(fNode);
        os.setFIndex(fIdx);
        for (int ft = 0; ft < nComp; ft++) {
            double fNorm = fNode.getComponentNorm(ft);
            if (fNorm < MachineZero) { continue; }
            os.setFComponent(ft);
            for (int gt = 0; gt < nComp; gt++) {
                os.setGComponent(gt);
                applyOperator(os);
            }
        }
    }
    // Coordinate scaling normalization in applied direction
    const double scaling_factor = std::pow(gNode.getMWTree().getMRA().getWorldBox().getScalingFactor(this->applyDir), oper->getOrder());
    for (int i = 0; i < gNode.getNCoefs(); i++) gNode.getCoefs()[i] /= scaling_factor;
    this->calc_t[mrcpp_get_thread_num()].stop();

    this->norm_t[mrcpp_get_thread_num()].resume();
    gNode.calcNorms();
    this->norm_t[mrcpp_get_thread_num()].stop();
}

/**
 * @brief Build the **operator band** of input nodes that influence a given output node.
 *
 * @param gNode   Output node for which to gather contributing input nodes.
 * @param idx_band Output list of input node indices (aligned with returned vector).
 * @return Vector of pointers to input nodes in @p fTree that lie within the
 *         operator bandwidth along `applyDir`.
 *
 * @details
 * The band extends `width = oper->getMaxBandWidth()` cells to the left/right
 * of the output index along the applied direction. Out-of-bounds indices are
 * skipped; periodicity is handled by FunctionTree::getRootIndex().
 */
template <int D, typename T> MWNodeVector<D, T> DerivativeCalculator<D, T>::makeOperBand(const MWNode<D, T> &gNode, std::vector<NodeIndex<D>> &idx_band) {
    assert(this->applyDir >= 0);
    assert(this->applyDir < D);

    MWNodeVector<D, T> band;
    const NodeIndex<D> &idx_0 = gNode.getNodeIndex();

    // Assumes given width only in applyDir, otherwise width = 0
    int width = this->oper->getMaxBandWidth();
    for (int w = -width; w <= width; w++) {
        NodeIndex<D> idx_w(idx_0);
        idx_w[this->applyDir] += w;

        // returns -1 if out of bounds and 0 for periodic
        int rIdx_w = this->fTree->getRootIndex(idx_w);
        if (rIdx_w >= 0) {
            idx_band.push_back(idx_w);
            band.push_back(&this->fTree->getNode(idx_w));
        }
    }
    return band;
}

/**
 * @brief Apply a single **zero-bandwidth** operator component to one input node.
 *
 * @param os Operator state (holds pointers/scratch and component indices).
 *
 * @details
 * Fetches the operator block at the current depth with translation 0 in all
 * directions. In non-applied directions, activates identity by passing
 * `nullptr` in the operator pointers, which signals tensorApplyOperComp() to
 * copy/accumulate.
 */
template <int D, typename T> void DerivativeCalculator<D, T>::applyOperator_bw0(OperatorState<D, T> &os) {
    MWNode<D, T> &gNode = *os.gNode;
    MWNode<D, T> &fNode = *os.fNode;
    const NodeIndex<D> &fIdx = *os.fIdx;
    const NodeIndex<D> &gIdx = gNode.getNodeIndex();
    int depth = gNode.getDepth();

    double oNorm = 1.0;
    double **oData = os.getOperData();

    for (int d = 0; d < D; d++) {
        const OperatorTree &oTree = this->oper->getComponent(0, d);
        const OperatorNode &oNode = oTree.getNode(depth, 0);
        int oIdx = os.getOperIndex(d);
        if (this->applyDir == d) {
            oData[d] = const_cast<double *>(oNode.getCoefs()) + oIdx * os.kp1_2;
        } else {
            if (oIdx == 0 or oIdx == 3) {
                // Identity in direction d
                oData[d] = nullptr;
            } else {
                // Outside identity block: contributes zero
                return;
            }
        }
    }
    this->operStat.incrementFNodeCounters(fNode, os.ft, os.gt);
    tensorApplyOperComp(os);
}

/**
 * @brief Apply a single **finite-bandwidth** operator component to one input node.
 *
 * @param os Operator state (holds pointers/scratch and component indices).
 *
 * @details
 * For each dimension:
 *  - Determine the relative translation from input to output node.
 *  - Check that translation lies within the operator bandwidth.
 *  - Fetch the corresponding #mrcpp::OperatorNode data.
 *  - In non-applied directions, only the central identity block (translation 0,
 *    component 0 or 3) contributes; otherwise the term is skipped.
 */
template <int D, typename T> void DerivativeCalculator<D, T>::applyOperator(OperatorState<D, T> &os) {
    MWNode<D, T> &gNode = *os.gNode;
    MWNode<D, T> &fNode = *os.fNode;
    const NodeIndex<D> &fIdx = *os.fIdx;
    const NodeIndex<D> &gIdx = gNode.getNodeIndex();
    int depth = gNode.getDepth();

    double oNorm = 1.0;
    double **oData = os.getOperData();

    for (int d = 0; d < D; d++) {
        const OperatorTree &oTree = this->oper->getComponent(0, d);

        int oTransl = fIdx[d] - gIdx[d];

        // Bandwidth check in each direction
        int a = (os.gt & (1 << d)) >> d;
        int b = (os.ft & (1 << d)) >> d;
        int idx = (a << 1) + b;
        int w = oTree.getBandWidth().getWidth(depth, idx);
        if (abs(oTransl) > w) { return; }

        const OperatorNode &oNode = oTree.getNode(depth, oTransl);
        int oIdx = os.getOperIndex(d);
        double ocn = oNode.getComponentNorm(oIdx);
        oNorm *= ocn;
        if (this->applyDir == d) {
            oData[d] = const_cast<double *>(oNode.getCoefs()) + oIdx * os.kp1_2;
        } else {
            if (oTransl == 0 and (oIdx == 0 or oIdx == 3)) {
                // Identity in direction d
                oData[d] = nullptr;
            } else {
                // Zero contribution
                return;
            }
        }
    }
    this->operStat.incrementFNodeCounters(fNode, os.ft, os.gt);
    tensorApplyOperComp(os);
}

/**
 * @brief Perform the separated **tensor contraction** for one operator term.
 *
 * @param os Operator state (provides temporary buffers and operator slices).
 *
 * @details
 * For each dimension i:
 *  - Map the \f$k\times k^{D-1}\f$ slice of the input into `f`,
 *  - Multiply by the \f$k\times k\f$ operator block if present
 *    (otherwise use identity),
 *  - Transpose-accumulate into the next staging buffer `g`.
 * On the last dimension, accumulate into the output buffer.
 */
template <int D, typename T> void DerivativeCalculator<D, T>::tensorApplyOperComp(OperatorState<D, T> &os) {
    T **aux = os.getAuxData();
    double **oData = os.getOperData();
    for (int i = 0; i < D; i++) {
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> f(aux[i], os.kp1, os.kp1_dm1);
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> g(aux[i + 1], os.kp1_dm1, os.kp1);
        if (oData[i] != nullptr) {
            Eigen::Map<MatrixXd> op(oData[i], os.kp1, os.kp1);
            if (i == D - 1) { // last dim: accumulate
                g.noalias() += f.transpose() * op;
            } else {
                g.noalias() = f.transpose() * op;
            }
        } else {
            // Identity in dimension i
            if (i == D - 1) {
                g.noalias() += f.transpose();
            } else {
                g.noalias() = f.transpose();
            }
        }
    }
}

/**
 * @brief Provide the initial work vector for a tree traversal.
 *
 * @param tree Output tree where results will be stored.
 * @return A vector of pointers to the leaf/end nodes of @p tree.
 *
 * @details
 * The derivative application uses a fixed grid determined by the operator.
 * This helper asks the tree to provide a snapshot of its end-node table to
 * seed the traversal.
 */
template <int D, typename T> MWNodeVector<D, T> *DerivativeCalculator<D, T>::getInitialWorkVector(MWTree<D, T> &tree) const {
    return tree.copyEndNodeTable();
}

// Explicit instantiations
template class DerivativeCalculator<1, double>;
template class DerivativeCalculator<2, double>;
template class DerivativeCalculator<3, double>;

template class DerivativeCalculator<1, ComplexDouble>;
template class DerivativeCalculator<2, ComplexDouble>;
template class DerivativeCalculator<3, ComplexDouble>;

} // namespace mrcpp
