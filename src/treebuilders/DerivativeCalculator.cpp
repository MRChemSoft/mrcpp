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

template <int D, typename T>
DerivativeCalculator<D, T>::DerivativeCalculator(int dir, DerivativeOperator<D> &o, FunctionTree<D, T> &f)
        : applyDir(dir)
        , fTree(&f)
        , oper(&o) {
    if (dir < 0 or dir >= D) MSG_ABORT("Invalid apply dir");
    initTimers();
}

template <int D, typename T> DerivativeCalculator<D, T>::~DerivativeCalculator() {
    this->operStat.flushNodeCounters();
    println(10, this->operStat);
}

template <int D, typename T> void DerivativeCalculator<D, T>::initTimers() {
    int nThreads = mrcpp_get_max_threads();
    for (int i = 0; i < nThreads; i++) {
        this->band_t.push_back(Timer(false));
        this->calc_t.push_back(Timer(false));
        this->norm_t.push_back(Timer(false));
    }
}

template <int D, typename T> void DerivativeCalculator<D, T>::clearTimers() {
    this->band_t.clear();
    this->calc_t.clear();
    this->norm_t.clear();
}

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

template <int D, typename T> void DerivativeCalculator<D, T>::calcNode(MWNode<D, T> &inpNode, MWNode<D, T> &outNode) {
    // if (this->oper->getMaxBandWidth() > 1) MSG_ABORT("Only implemented for zero bw");
    outNode.zeroCoefs();
    int nComp = (1 << D);
    std::vector<T> tmpCoefs(outNode.getNCoefs());
    OperatorState<D, T> os(outNode, tmpCoefs.data());

    os.setFNode(inpNode);
    os.setFIndex(inpNode.nodeIndex);
    for (int ft = 0; ft < nComp; ft++) {
        double fNorm = inpNode.getComponentNorm(ft);
        if (fNorm < MachineZero) { continue; } // Could this be used outside the loop?
        os.setFComponent(ft);
        for (int gt = 0; gt < nComp; gt++) {
            os.setGComponent(gt);
            applyOperator_bw0(os);
        }
    }
    // Multiply appropriate scaling factor. TODO: Could be included elsewhere
    const double scaling_factor = 1.0 / std::pow(outNode.getMWTree().getMRA().getWorldBox().getScalingFactor(this->applyDir), oper->getOrder());
    if (abs(scaling_factor - 1.0) > MachineZero) {
        for (int i = 0; i < outNode.getNCoefs(); i++) outNode.getCoefs()[i] *= scaling_factor;
    }
    outNode.calcNorms(); // TODO:required? norms are not used for now
}

template <int D, typename T> void DerivativeCalculator<D, T>::calcNode(MWNode<D, T> &gNode) {
    gNode.zeroCoefs();

    int nComp = (1 << D);
    std::vector<T> tmpCoefs(gNode.getNCoefs());
    OperatorState<D, T> os(gNode, tmpCoefs.data());
    this->operStat.incrementGNodeCounters(gNode);

    // Get all nodes in f within the bandwith of O in g
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
    // Multiply appropriate scaling factor
    const double scaling_factor = std::pow(gNode.getMWTree().getMRA().getWorldBox().getScalingFactor(this->applyDir), oper->getOrder());
    for (int i = 0; i < gNode.getNCoefs(); i++) gNode.getCoefs()[i] /= scaling_factor;
    this->calc_t[mrcpp_get_thread_num()].stop();

    this->norm_t[mrcpp_get_thread_num()].resume();
    gNode.calcNorms();
    this->norm_t[mrcpp_get_thread_num()].stop();
}

/** Return a vector of nodes in F affected by O, given a node in G */
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

/** Apply a single operator component (term) to a single f-node assuming zero bandwidth */
template <int D, typename T> void DerivativeCalculator<D, T>::applyOperator_bw0(OperatorState<D, T> &os) {
    // cout<<" applyOperator "<<endl;
    MWNode<D, T> &gNode = *os.gNode;
    MWNode<D, T> &fNode = *os.fNode;
    int depth = gNode.getDepth();

    double **oData = os.getOperData();

    for (int d = 0; d < D; d++) {
        const OperatorTree &oTree = this->oper->getComponent(0, d);
        const OperatorNode &oNode = oTree.getNode(depth, 0);
        int oIdx = os.getOperIndex(d);
        if (this->applyDir == d) {
            oData[d] = const_cast<double *>(oNode.getCoefs()) + oIdx * os.kp1_2;
        } else {
            if (oIdx == 0 or oIdx == 3) {
                // This will activate the identity operator in direction i
                oData[d] = nullptr;
            } else {
                // This means that we are in a zero part of the identity operator
                return;
            }
        }
    }
    this->operStat.incrementFNodeCounters(fNode, os.ft, os.gt);
    tensorApplyOperComp(os);
}

/** Apply a single operator component (term) to a single f-node. Whether the
operator actualy is applied is determined by a screening threshold. */
template <int D, typename T> void DerivativeCalculator<D, T>::applyOperator(OperatorState<D, T> &os) {
    MWNode<D, T> &gNode = *os.gNode;
    MWNode<D, T> &fNode = *os.fNode;
    const NodeIndex<D> &fIdx = *os.fIdx;
    const NodeIndex<D> &gIdx = gNode.getNodeIndex();
    int depth = gNode.getDepth();

    double **oData = os.getOperData();

    for (int d = 0; d < D; d++) {
        const OperatorTree &oTree = this->oper->getComponent(0, d);

        int oTransl = fIdx[d] - gIdx[d];

        //  The following will check the actual band width in each direction.
        //  Not needed if the thresholding at the end of this routine is active.
        int a = (os.gt & (1 << d)) >> d;
        int b = (os.ft & (1 << d)) >> d;
        int idx = (a << 1) + b;
        int w = oTree.getBandWidth().getWidth(depth, idx);
        if (abs(oTransl) > w) { return; }

        const OperatorNode &oNode = oTree.getNode(depth, oTransl);
        int oIdx = os.getOperIndex(d);
        if (this->applyDir == d) {
            oData[d] = const_cast<double *>(oNode.getCoefs()) + oIdx * os.kp1_2;
        } else {
            if (oTransl == 0 and (oIdx == 0 or oIdx == 3)) {
                // This will activate the identity operator in direction i
                oData[d] = nullptr;
            } else {
                // This means that we are in a zero part of the identity operator
                return;
            }
        }
    }
    this->operStat.incrementFNodeCounters(fNode, os.ft, os.gt);
    tensorApplyOperComp(os);
}

/** Perform the required linear algebra operations in order to apply an
operator component to a f-node in a n-dimensional tensor space. */
template <int D, typename T> void DerivativeCalculator<D, T>::tensorApplyOperComp(OperatorState<D, T> &os) {
    T **aux = os.getAuxData();
    double **oData = os.getOperData();
    for (int i = 0; i < D; i++) {
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> f(aux[i], os.kp1, os.kp1_dm1);
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> g(aux[i + 1], os.kp1_dm1, os.kp1);
        if (oData[i] != nullptr) {
            Eigen::Map<MatrixXd> op(oData[i], os.kp1, os.kp1);
            if (i == D - 1) { // Last dir: Add up into g
                g.noalias() += f.transpose() * op;
            } else {
                g.noalias() = f.transpose() * op;
            }
        } else {
            // Identity operator in direction i
            if (i == D - 1) { // Last dir: Add up into g
                g.noalias() += f.transpose();
            } else {
                g.noalias() = f.transpose();
            }
        }
    }
}

template <int D, typename T> MWNodeVector<D, T> *DerivativeCalculator<D, T>::getInitialWorkVector(MWTree<D, T> &tree) const {
    return tree.copyEndNodeTable();
}

template class DerivativeCalculator<1, double>;
template class DerivativeCalculator<2, double>;
template class DerivativeCalculator<3, double>;

template class DerivativeCalculator<1, ComplexDouble>;
template class DerivativeCalculator<2, ComplexDouble>;
template class DerivativeCalculator<3, ComplexDouble>;

} // namespace mrcpp
