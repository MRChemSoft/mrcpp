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

#include "ConvolutionCalculator.h"
#include "operators/ConvolutionOperator.h"
#include "operators/OperatorState.h"
#include "trees/BandWidth.h"
#include "trees/FunctionNode.h"
#include "trees/OperatorNode.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/tree_utils.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using Eigen::MatrixXd;
using Eigen::MatrixXi;

namespace mrcpp {

template <int D>
ConvolutionCalculator<D>::ConvolutionCalculator(double p, ConvolutionOperator<D> &o, FunctionTree<D> &f, int depth)
        : maxDepth(depth)
        , prec(p)
        , oper(&o)
        , fTree(&f) {
    if (this->maxDepth > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    initBandSizes();
    initTimers();
}

template <int D> ConvolutionCalculator<D>::~ConvolutionCalculator() {
    clearTimers();
    this->operStat.flushNodeCounters();
    println(10, this->operStat);
    for (int i = 0; i < this->bandSizes.size(); i++) { delete this->bandSizes[i]; }
}

template <int D> void ConvolutionCalculator<D>::initTimers() {
    int nThreads = mrcpp_get_max_threads();
    for (int i = 0; i < nThreads; i++) {
        this->band_t.push_back(new Timer(false));
        this->calc_t.push_back(new Timer(false));
        this->norm_t.push_back(new Timer(false));
    }
}

template <int D> void ConvolutionCalculator<D>::clearTimers() {
    int nThreads = mrcpp_get_max_threads();
    for (int i = 0; i < nThreads; i++) {
        delete this->band_t[i];
        delete this->calc_t[i];
        delete this->norm_t[i];
    }
    this->band_t.clear();
    this->calc_t.clear();
    this->norm_t.clear();
}

template <int D> void ConvolutionCalculator<D>::printTimers() const {
    int oldprec = Printer::setPrecision(1);
    int nThreads = mrcpp_get_max_threads();
    printout(20, "\n\nthread ");
    for (int i = 0; i < nThreads; i++) printout(20, std::setw(9) << i);
    printout(20, "\nband     ");
    for (int i = 0; i < nThreads; i++) printout(20, this->band_t[i]->elapsed() << "  ");
    printout(20, "\ncalc     ");
    for (int i = 0; i < nThreads; i++) printout(20, this->calc_t[i]->elapsed() << "  ");
    printout(20, "\nnorm     ");
    for (int i = 0; i < nThreads; i++) printout(20, this->norm_t[i]->elapsed() << "  ");
    printout(20, "\n\n");
    Printer::setPrecision(oldprec);
}

/** Initialize the number of nodes formally within the bandwidth of an
 operator. The band size is used for thresholding. */
template <int D> void ConvolutionCalculator<D>::initBandSizes() {
    for (int i = 0; i < this->oper->size(); i++) {
        const OperatorTree &oTree = *(*this->oper)[i];
        const BandWidth &bw = oTree.getBandWidth();
        auto *bsize = new MatrixXi(this->maxDepth, this->nComp2 + 1);
        bsize->setZero();
        for (int j = 0; j < this->maxDepth; j++) { calcBandSizeFactor(*bsize, j, bw); }
        this->bandSizes.push_back(bsize);
    }
}

/** Calculate the number of nodes within the bandwidth
 * of an operator. Currently this routine ignores the fact that
 * there are edges on the world box, and thus over estimates
 * the number of nodes. This is different from the previous version. */
template <int D> void ConvolutionCalculator<D>::calcBandSizeFactor(MatrixXi &bs, int depth, const BandWidth &bw) {
    for (int gt = 0; gt < this->nComp; gt++) {
        for (int ft = 0; ft < this->nComp; ft++) {
            int k = gt * this->nComp + ft;
            int totNodes = 1;
            for (int i = 0; i < D; i++) {
                int oIdx = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                int width = bw.getWidth(depth, oIdx);
                if (width < 0) {
                    bs(depth, k) = 0;
                    continue;
                }
                totNodes *= 2 * width + 1;
            }
            bs(depth, k) = totNodes * this->nComp2;
        }
    }
    bs(depth, this->nComp2) = bs.row(depth).maxCoeff();
}

/** Return a vector of nodes in F affected by O, given a node in G */
template <int D>
MWNodeVector<D> *ConvolutionCalculator<D>::makeOperBand(const MWNode<D> &gNode, std::vector<NodeIndex<D>> &idx_band) {
    auto *band = new MWNodeVector<D>;

    int depth = gNode.getDepth();
    bool periodic = gNode.getMWTree().getMRA().getWorldBox().isPeriodic();
    int width = this->oper->getMaxBandWidth(depth);
    if (width >= 0) {
        const NodeBox<D> &fWorld = this->fTree->getRootBox();
        const NodeIndex<D> &cIdx = fWorld.getCornerIndex();
        const NodeIndex<D> &gIdx = gNode.getNodeIndex();

        int nbox[D];
        int scale = gNode.getScale();
        NodeIndex<D> sIdx(scale); // start index
        NodeIndex<D> eIdx(scale); // end index
        for (int i = 0; i < D; i++) {
            sIdx[i] = gIdx[i] - width;
            eIdx[i] = gIdx[i] + width;
            // We need to consider the world borders
            int nboxes = fWorld.size(i) * (1 << depth);
            int c_i = cIdx[i] * (1 << depth);
            if (sIdx[i] < c_i and !periodic) sIdx[i] = c_i;
            if (eIdx[i] > c_i + nboxes - 1 and !periodic) eIdx[i] = c_i + nboxes - 1;
            nbox[i] = eIdx[i] - sIdx[i] + 1;
        }

        fillOperBand(band, idx_band, sIdx, nbox, D - 1);
    }
    return band;
}

/** Recursively retrieve all reachable f-nodes within the bandwidth. */
template <int D>
void ConvolutionCalculator<D>::fillOperBand(MWNodeVector<D> *band,
                                            std::vector<NodeIndex<D>> &idx_band,
                                            NodeIndex<D> &idx,
                                            const int *nbox,
                                            int dim) {
    int l_start = idx[dim];
    for (int j = 0; j < nbox[dim]; j++) {
        // Recurse until dim == 0
        if (dim > 0) {
            fillOperBand(band, idx_band, idx, nbox, dim - 1);
            idx[dim]++;
            continue;
        }
        MWNode<D> &fNode = this->fTree->getNode(idx);
        idx_band.push_back(idx);
        band->push_back(&fNode);
        idx[dim]++;
    }
    idx[dim] = l_start;
}

template <int D> void ConvolutionCalculator<D>::calcNode(MWNode<D> &node) {
    auto &gNode = static_cast<FunctionNode<D> &>(node);
    gNode.zeroCoefs();

    int depth = gNode.getDepth();
    double tmpCoefs[gNode.getNCoefs()];
    OperatorState<D> os(gNode, tmpCoefs);
    this->operStat.incrementGNodeCounters(gNode);

    // Get all nodes in f within the bandwith of O in g
    this->band_t[mrcpp_get_thread_num()]->resume();
    std::vector<NodeIndex<D>> idx_band;
    MWNodeVector<D> *fBand = makeOperBand(gNode, idx_band);
    this->band_t[mrcpp_get_thread_num()]->stop();

    MWTree<D> &gTree = gNode.getMWTree();
    double gThrs = gTree.getSquareNorm();
    if (gThrs > 0.0) {
        auto nTerms = static_cast<double>(this->oper->size());
        auto precFac = this->precFunc(gNode.getNodeIndex());
        gThrs = this->prec * precFac * std::sqrt(gThrs / nTerms);
    }

    os.gThreshold = gThrs;

    this->calc_t[mrcpp_get_thread_num()]->resume();
    for (int n = 0; n < fBand->size(); n++) {
        MWNode<D> &fNode = *(*fBand)[n];
        NodeIndex<D> &fIdx = idx_band[n];
        os.setFNode(fNode);
        os.setFIndex(fIdx);
        for (int ft = 0; ft < this->nComp; ft++) {
            double fNorm = fNode.getComponentNorm(ft);
            if (fNorm < MachineZero) { continue; }
            os.setFComponent(ft);
            for (int gt = 0; gt < this->nComp; gt++) {
                if (depth == 0 or gt != 0 or ft != 0) {
                    os.setGComponent(gt);
                    applyOperComp(os);
                }
            }
        }
    }
    this->calc_t[mrcpp_get_thread_num()]->stop();

    this->norm_t[mrcpp_get_thread_num()]->resume();
    gNode.calcNorms();
    this->norm_t[mrcpp_get_thread_num()]->stop();
    delete fBand;
}

/** Apply each component (term) of the operator expansion to a node in f */
template <int D> void ConvolutionCalculator<D>::applyOperComp(OperatorState<D> &os) {
    int depth = os.gNode->getDepth();
    double fNorm = os.fNode->getComponentNorm(os.ft);
    for (int i = 0; i < this->oper->size(); i++) {
        const OperatorTree &ot = *(*this->oper)[i];
        const BandWidth &bw = ot.getBandWidth();
        if (os.getMaxDeltaL() > bw.getMaxWidth(depth)) { continue; }
        os.oTree = &ot;
        os.fThreshold = getBandSizeFactor(i, depth, os) * fNorm;
        applyOperator(os);
    }
}

/** Apply a single operator component (term) to a single f-node. Whether the
operator actualy is applied is determined by a screening threshold. */
template <int D> void ConvolutionCalculator<D>::applyOperator(OperatorState<D> &os) {
    const OperatorTree &oTree = *os.oTree;
    MWNode<D> &gNode = *os.gNode;
    MWNode<D> &fNode = *os.fNode;
    const NodeIndex<D> &fIdx = *os.fIdx;
    const NodeIndex<D> &gIdx = gNode.getNodeIndex();
    int depth = gNode.getDepth();

    double oNorm = 1.0;
    double **oData = os.getOperData();

    for (int d = 0; d < D; d++) {
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
        oNorm *= oNode.getComponentNorm(oIdx);
        oData[d] = const_cast<double *>(oNode.getCoefs()) + oIdx * os.kp1_2;
    }
    double upperBound = oNorm * os.fThreshold;
    if (upperBound > os.gThreshold) {
        this->operStat.incrementFNodeCounters(fNode, os.ft, os.gt);
        tensorApplyOperComp(os);
    }
}

/** Perorm the required linear algebra operations in order to apply an
operator component to a f-node in a n-dimensional tesor space. */
template <int D> void ConvolutionCalculator<D>::tensorApplyOperComp(OperatorState<D> &os) {
    double **aux = os.getAuxData();
    double **oData = os.getOperData();
#ifdef HAVE_BLAS
    double mult = 0.0;
    for (int i = 0; i < D; i++) {
        if (oData[i] != 0) {
            if (i == D - 1) { // Last dir: Add up into g
                mult = 1.0;
            }
            const double *f = aux[i];
            double *g = const_cast<double *>(aux[i + 1]);
            cblas_dgemm(CblasColMajor,
                        CblasTrans,
                        CblasNoTrans,
                        os.kp1_dm1,
                        os.kp1,
                        os.kp1,
                        1.0,
                        f,
                        os.kp1,
                        oData[i],
                        os.kp1,
                        mult,
                        g,
                        os.kp1_dm1);
        } else {
            // Identity operator in direction i
            Eigen::Map<MatrixXd> f(aux[i], os.kp1, os.kp1_dm1);
            Eigen::Map<MatrixXd> g(aux[i + 1], os.kp1_dm1, os.kp1);
            if (oData[i] == 0) {
                if (i == D - 1) { // Last dir: Add up into g
                    g += f.transpose();
                } else {
                    g = f.transpose();
                }
            }
        }
    }
#else
    for (int i = 0; i < D; i++) {
        Eigen::Map<MatrixXd> f(aux[i], os.kp1, os.kp1_dm1);
        Eigen::Map<MatrixXd> g(aux[i + 1], os.kp1_dm1, os.kp1);
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
#endif
}

template <int D> MWNodeVector<D> *ConvolutionCalculator<D>::getInitialWorkVector(MWTree<D> &tree) const {
    auto *nodeVec = new MWNodeVector<D>;
    tree_utils::make_node_table(tree, *nodeVec);
    return nodeVec;
}

template class ConvolutionCalculator<1>;
template class ConvolutionCalculator<2>;
template class ConvolutionCalculator<3>;

} // namespace mrcpp
