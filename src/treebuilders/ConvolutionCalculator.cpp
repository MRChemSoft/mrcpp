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
 * @file ConvolutionCalculator.cpp
 * @brief Adaptive node-wise application kernel for separable convolution operators.
 *
 * @details
 * This file implements the templated class
 * mrcpp::ConvolutionCalculator, which is the **workhorse** used by the
 * adaptive `TreeBuilder` when applying a separable convolution operator
 * (#mrcpp::ConvolutionOperator) to a multiresolution function tree
 * (#mrcpp::FunctionTree).
 *
 * At a high level, for each **target** node \f$ g \f$ (in the output tree)
 * the calculator:
 *  - determines the **band** of **source** nodes \f$ f \f$ that can
 *    contribute via the operator's bandwidth model,
 *  - estimates cheap **screening bounds** using precomputed operator norms,
 *    the local source/target norms, and a precision policy,
 *  - for surviving pairs \f$ (g,f) \f$, performs a sequence of small
 *    **tensor contractions** (one per Cartesian direction) to apply the
 *    separable operator component(s) and accumulates the result into \f$ g \f$.
 *
 * The class also:
 *  - precomputes **band-size factors** per depth and component-combination to
 *    drive thresholding,
 *  - supports **periodic worlds** and optional **unit-cell manipulation**
 *    (near-field vs. far-field selection),
 *  - collects **per-thread timings** and **operator-usage statistics**.
 *
 * ### Screening model (outline)
 * Let \f$ \mathcal{O} = \sum_i \bigotimes_{d=1}^D O_i^{(d)} \f$ be the
 * separable expansion (terms indexed by \f$ i \f$). For a source node
 * \f$ f \f$ and target node \f$ g \f$, the calculator estimates
 * \f[
 *   \| \mathcal{O}_i f \| \;\lesssim\;
 *   \Big(\prod_{d=1}^D \|O_i^{(d)}\|\Big)\; \|f\|\; s(i, \Delta \ell)
 * \f]
 * where \f$ s(\cdot) \f$ is a band-size factor depending on depth and the
 * component combination, and compares the bound to a target threshold
 * \f$ \tau(g) \sim \texttt{prec} \cdot \sqrt{\|g\|^2 / N_\text{terms}} \f$.
 * Only terms that can exceed \f$ \tau(g) \f$ are explicitly applied.
 *
 * ### BLAS vs. Eigen
 * If BLAS is available, the directional contractions can be carried out via
 * GEMM. Otherwise, an Eigen-based path is used. Both routes compute
 * \f$ G \leftarrow F^\top O \f$ in each direction and **accumulate** on the
 * last direction to the target buffer.
 */

#include "ConvolutionCalculator.h"
#include "operators/ConvolutionOperator.h"
#include "operators/OperatorState.h"
#include "trees/BandWidth.h"
#include "trees/FunctionNode.h"
#include "trees/OperatorNode.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"
#include "utils/periodic_utils.h"
#include "utils/tree_utils.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using Eigen::MatrixXd;
using Eigen::MatrixXi;

namespace mrcpp {

/**
 * @brief Construct a calculator for applying a convolution operator.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (`double` or `ComplexDouble`).
 * @param p       Target precision used for screening and adaptivity.
 * @param o       Separable convolution operator to apply.
 * @param f       Source function tree (input).
 * @param depth   Maximum operator depth considered for band-size tables.
 *
 * @details
 * Initializes per-term **band-size tables** (used in screening) and
 * allocates per-thread timers. The `depth` argument is upper-bounded by
 * `MaxDepth`.
 */
template <int D, typename T>
ConvolutionCalculator<D, T>::ConvolutionCalculator(double p, ConvolutionOperator<D> &o, FunctionTree<D, T> &f, int depth)
        : maxDepth(depth)
        , prec(p)
        , oper(&o)
        , fTree(&f) {
    if (this->maxDepth > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    initBandSizes();
    initTimers();
}

/**
 * @brief Destructor: clear timers and print aggregated operator statistics.
 */
template <int D, typename T> ConvolutionCalculator<D, T>::~ConvolutionCalculator() {
    clearTimers();
    this->operStat.flushNodeCounters();
    println(10, this->operStat);
    for (int i = 0; i < this->bandSizes.size(); i++) { delete this->bandSizes[i]; }
}

/**
 * @brief Allocate per-thread timers for band construction, calculation, and norm updates.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::initTimers() {
    int nThreads = mrcpp_get_max_threads();
    for (int i = 0; i < nThreads; i++) {
        this->band_t.push_back(new Timer(false));
        this->calc_t.push_back(new Timer(false));
        this->norm_t.push_back(new Timer(false));
    }
}

/**
 * @brief Release per-thread timers.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::clearTimers() {
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

/**
 * @brief Print a compact report of thread-wise timings.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::printTimers() const {
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

/**
 * @brief Precompute per-depth band-size factors for all operator terms.
 *
 * @details
 * For each raw operator term and each depth, builds a table of the number of
 * source nodes formally falling within the **Cartesian bandwidth box** for
 * every component-combination (gt,ft). These factors are later used to scale
 * screening thresholds.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::initBandSizes() {
    for (int i = 0; i < this->oper->size(); i++) {
        // IMPORTANT: only 0-th dimension!
        const OperatorTree &oTree = this->oper->getComponent(i, 0);
        const BandWidth &bw = oTree.getBandWidth();
        auto *bsize = new MatrixXi(this->maxDepth, this->nComp2 + 1);
        bsize->setZero();
        for (int j = 0; j < this->maxDepth; j++) { calcBandSizeFactor(*bsize, j, bw); }
        this->bandSizes.push_back(bsize);
    }
}

/**
 * @brief Compute band-size factor for a given depth from a bandwidth model.
 *
 * @param[out] bs   Table to be filled (rows: depth, cols: component-pairs plus a max column).
 * @param[in]  depth Operator depth relative to root.
 * @param[in]  bw   Bandwidth model (per-depth widths per component index).
 *
 * @details
 * For each component pair \f$(g_t,f_t)\f$, the routine forms the Cartesian
 * product of directional half-widths to estimate the number of contributing
 * source nodes and stores it in \p bs. The last column stores the row-wise
 * maximum for quick access.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::calcBandSizeFactor(MatrixXi &bs, int depth, const BandWidth &bw) {
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

/**
 * @brief Build the band of source nodes affected by the operator for a given target node.
 *
 * @param[in]  gNode    Target node (in the output tree).
 * @param[out] idx_band Matching indices of the source nodes added to the band.
 * @returns A vector of pointers to the source nodes \f$ f \f$.
 *
 * @details
 * The band is the intersection between the operator's bandwidth box centered
 * at \p gNode and the function-tree world box, respecting periodicity and
 * (optionally) unit-cell filtering when `manipulateOperator` is enabled.
 */
template <int D, typename T> MWNodeVector<D, T> *ConvolutionCalculator<D, T>::makeOperBand(const MWNode<D, T> &gNode, std::vector<NodeIndex<D>> &idx_band) {
    auto *band = new MWNodeVector<D, T>;

    int o_depth = gNode.getScale() - this->oper->getOperatorRoot();
    int g_depth = gNode.getDepth();
    int width = this->oper->getMaxBandWidth(o_depth);

    bool periodic = gNode.getMWTree().isPeriodic();
    int reach = this->oper->getOperatorReach();

    if (width >= 0) {
        const NodeBox<D, T> &fWorld = this->fTree->getRootBox();
        const NodeIndex<D> &cIdx = fWorld.getCornerIndex();
        const NodeIndex<D> &gIdx = gNode.getNodeIndex();

        int nbox[D];
        int scale = gNode.getScale();
        NodeIndex<D> sIdx(scale); // start index
        NodeIndex<D> eIdx(scale); // end index
        for (int i = 0; i < D; i++) {
            sIdx[i] = gIdx[i] - width;
            eIdx[i] = gIdx[i] + width;
            // Consider world borders / periodic wrapping
            int nboxes = fWorld.size(i) * (1 << o_depth);
            int c_i = cIdx[i] * (1 << o_depth);
            if (not periodic) {
                if (sIdx[i] < c_i) sIdx[i] = c_i;
                if (eIdx[i] > c_i + nboxes - 1) eIdx[i] = c_i + nboxes - 1;
            } else {
                if (sIdx[i] < c_i * reach) sIdx[i] = c_i * reach;
                if (eIdx[i] > (c_i + nboxes) * reach - 1) eIdx[i] = (c_i + nboxes) * reach - 1;
            }
            nbox[i] = eIdx[i] - sIdx[i] + 1;
        }

        fillOperBand(band, idx_band, sIdx, nbox, D - 1);
    }
    return band;
}

/**
 * @brief Recursive helper to enumerate all source indices inside the bandwidth box.
 *
 * @param[out] band     Vector of pointers to source nodes added along the recursion.
 * @param[out] idx_band Parallel vector of node indices corresponding to \p band.
 * @param[in]  idx      Current multi-index (mutated along recursion).
 * @param[in]  nbox     Side lengths of the bandwidth box.
 * @param[in]  dim      Current dimension to recurse on.
 *
 * @details
 * If **unit-cell manipulation** is enabled, nodes are included/excluded based
 * on their membership in the first unit cell (for periodic worlds) and the
 * `onUnitcell` flag.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::fillOperBand(MWNodeVector<D, T> *band, std::vector<NodeIndex<D>> &idx_band, NodeIndex<D> &idx, const int *nbox, int dim) {
    int l_start = idx[dim];
    for (int j = 0; j < nbox[dim]; j++) {
        // Recurse until dim == 0
        if (dim > 0) {
            fillOperBand(band, idx_band, idx, nbox, dim - 1);
            idx[dim]++;
            continue;
        }
        if (not manipulateOperator) {
            MWNode<D, T> &fNode = this->fTree->getNode(idx);
            idx_band.push_back(idx);
            band->push_back(&fNode);

        } else {
            const auto oper_scale = this->oper->getOperatorRoot();
            if (oper_scale == 0) {
                if (periodic::in_unit_cell<D>(idx) and onUnitcell) {
                    MWNode<D, T> &fNode = this->fTree->getNode(idx);
                    idx_band.push_back(idx);
                    band->push_back(&fNode);
                }
                if (not periodic::in_unit_cell<D>(idx) and not onUnitcell) {
                    MWNode<D, T> &fNode = this->fTree->getNode(idx);
                    idx_band.push_back(idx);
                    band->push_back(&fNode);
                }
            } else if (oper_scale < 0) {
                if (periodic::in_unit_cell<D>(idx) and onUnitcell) {
                    MWNode<D, T> &fNode = this->fTree->getNode(idx);
                    idx_band.push_back(idx);
                    band->push_back(&fNode);
                }
                if (not onUnitcell) MSG_ABORT("Cannot do with negative operator scale");
            } else
                MSG_ABORT("Cannot manipulate operators with positive operator scale");
        }
        idx[dim]++;
    }
    idx[dim] = l_start;
}

/**
 * @brief Compute contributions to a single **target** node by scanning its band.
 *
 * @param[in,out] node Target node (coefficients are accumulated here).
 *
 * @details
 * - Builds the source band for the target node.
 * - Computes a **local target threshold** from the node's tree norm and `prec`.
 * - Loops over band nodes and component combinations, performing **screening**.
 * - For surviving pairs, applies all operator terms via `applyOperComp`.
 * - Updates node norms at the end.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::calcNode(MWNode<D, T> &node) {
    auto &gNode = static_cast<FunctionNode<D, T> &>(node);
    gNode.zeroCoefs();

    int o_depth = gNode.getScale() - this->oper->getOperatorRoot();
    if (manipulateOperator and this->oper->getOperatorRoot() < 0) o_depth = gNode.getDepth();
    T tmpCoefs[gNode.getNCoefs()];
    OperatorState<D, T> os(gNode, tmpCoefs);
    this->operStat.incrementGNodeCounters(gNode);

    // Get all nodes in f within the bandwidth of O around g
    this->band_t[mrcpp_get_thread_num()]->resume();
    std::vector<NodeIndex<D>> idx_band;
    MWNodeVector<D, T> *fBand = makeOperBand(gNode, idx_band);
    this->band_t[mrcpp_get_thread_num()]->stop();

    // Build target threshold (relative by default; may be scaled by precFunc)
    MWTree<D, T> &gTree = gNode.getMWTree();
    double gThrs = gTree.getSquareNorm();
    if (gThrs > 0.0) {
        auto nTerms = static_cast<double>(this->oper->size());
        auto precFac = this->precFunc(gNode.getNodeIndex());
        gThrs = this->prec * precFac * std::sqrt(gThrs / nTerms);
    }
    os.gThreshold = gThrs;

    // Scan band and apply screened operator terms
    this->calc_t[mrcpp_get_thread_num()]->resume();
    for (int n = 0; n < fBand->size(); n++) {
        MWNode<D, T> &fNode = *(*fBand)[n];
        NodeIndex<D> &fIdx = idx_band[n];
        os.setFNode(fNode);
        os.setFIndex(fIdx);
        for (int ft = 0; ft < this->nComp; ft++) {
            double fNorm = fNode.getComponentNorm(ft);
            if (fNorm < MachineZero) { continue; }
            os.setFComponent(ft);
            for (int gt = 0; gt < this->nComp; gt++) {
                if (o_depth == 0 or gt != 0 or ft != 0) {
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

/**
 * @brief Apply all operator **terms** for a fixed component pair (ft,gt), with screening.
 *
 * @param[in,out] os Transient operator state bundling node pointers, buffers and norms.
 *
 * @details
 * For each term \f$ i \f$:
 *  - Check per-depth **bandwidth** feasibility (via `BandWidth`).
 *  - Build a per-term **source threshold** from band-size tables and \f$ \|f\| \f$.
 *  - Compute an **upper bound** using cached directional norms of the operator
 *    node at the appropriate translation; compare to \f$ \tau(g) \f$.
 *  - If the bound exceeds \f$ \tau(g) \f$, execute the tensor contraction
 *    (see `tensorApplyOperComp`).
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::applyOperComp(OperatorState<D, T> &os) {
    double fNorm = os.fNode->getComponentNorm(os.ft);
    int o_depth = os.fNode->getScale() - this->oper->getOperatorRoot();
    for (int i = 0; i < this->oper->size(); i++) {
        // IMPORTANT: only 0-th dimension
        const OperatorTree &ot = this->oper->getComponent(i, 0);
        const BandWidth &bw = ot.getBandWidth();
        if (os.getMaxDeltaL() > bw.getMaxWidth(o_depth)) { continue; }
        os.fThreshold = getBandSizeFactor(i, o_depth, os) * fNorm;
        applyOperator(i, os);
    }
}

/**
 * @brief Apply a single operator term to a single source node (low-level path).
 *
 * @param i  Index of the operator term in the separable expansion.
 * @param os Operator state (nodes, buffers, norms, component indices).
 *
 * @details
 * For each direction:
 *  - Fetch the operator-block at the required translation (\f$ \Delta \ell \f$)
 *    and depth \f$ o\_depth \f$; multiply the running contraction with its norm
 *    and keep a raw pointer to its coefficient block.
 *  - If the translation is outside bandwidth, return early.
 * After the per-direction setup:
 *  - Form an **upper bound** as product of directional norms times the
 *    source-threshold and compare to the target-threshold.
 *  - If active, dispatch to `tensorApplyOperComp` to carry out the contraction
 *    and accumulate into the target node buffer.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::applyOperator(int i, OperatorState<D, T> &os) {
    MWNode<D, T> &gNode = *os.gNode;
    MWNode<D, T> &fNode = *os.fNode;

    const NodeIndex<D> &fIdx = *os.fIdx;
    const NodeIndex<D> &gIdx = gNode.getNodeIndex();
    int o_depth = gNode.getScale() - this->oper->getOperatorRoot();

    double oNorm = 1.0;
    double **oData = os.getOperData();

    for (int d = 0; d < D; d++) {
        auto &oTree = this->oper->getComponent(i, d);
        int oTransl = fIdx[d] - gIdx[d];

        // Per-direction bandwidth check
        int a = (os.gt & (1 << d)) >> d;
        int b = (os.ft & (1 << d)) >> d;
        int idx = (a << 1) + b;
        if (oTree.isOutsideBand(oTransl, o_depth, idx)) { return; }

        const OperatorNode &oNode = oTree.getNode(o_depth, oTransl);
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

/**
 * @brief Perform the directional tensor contractions for one operator term.
 *
 * @param os Operator state (holds mapped buffers for in-place contractions).
 *
 * @details
 * The contraction sequence computes, for each direction \f$ d \f$,
 * \f$ G \leftarrow F^\top O^{(d)} \f$, with **accumulation** on the last
 * direction. If a directional block is `nullptr`, an identity map is used
 * (i.e., pure transposition).
 *
 * Both a BLAS path (disabled here) and an Eigen path are implemented.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::tensorApplyOperComp(OperatorState<D, T> &os) {
    T **aux = os.getAuxData();
    double **oData = os.getOperData();
    /*
#ifdef HAVE_BLAS
    double mult = 0.0;
    for (int i = 0; i < D; i++) {
        if (oData[i] != 0) {
            if (i == D - 1) { // Last dir: Add up into g
                mult = 1.0;
            }
            const double *f = aux[i];
            double *g = const_cast<double *>(aux[i + 1]);
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, os.kp1_dm1, os.kp1, os.kp1, 1.0, f, os.kp1, oData[i], os.kp1, mult, g, os.kp1_dm1);
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
    */
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
    //#endif
}

/**
 * @brief Ensure parent nodes exist up to the operator root (periodic worlds).
 *
 * @param tree Target/output tree.
 *
 * @details
 * When operating in periodic settings, parent nodes above the root scale
 * may be required for coarse contributions; this helper guarantees their
 * presence prior to work scheduling.
 */
template <int D, typename T> void ConvolutionCalculator<D, T>::touchParentNodes(MWTree<D, T> &tree) const {
    if (not manipulateOperator) {
        const auto oper_scale = this->oper->getOperatorRoot();
        auto car_prod = math_utils::cartesian_product(std::vector<int>{-1, 0}, D);
        for (auto i = -1; i > oper_scale - 1; i--) {
            for (auto &a : car_prod) {
                std::array<int, D> l{};
                std::copy_n(a.begin(), D, l.begin());
                NodeIndex<D> idx(i, l);
                tree.getNode(idx);
                this->fTree->getNode(idx);
            }
        }
    }
}

/**
 * @brief Create the initial list of target nodes to process.
 *
 * @param tree Target/output tree.
 * @returns A vector of pointers to existing nodes to be processed.
 *
 * @details
 * For periodic trees, parent nodes above the root are first touched to ensure
 * consistency; then a flat node table is produced via `tree_utils::make_node_table`.
 */
template <int D, typename T> MWNodeVector<D, T> *ConvolutionCalculator<D, T>::getInitialWorkVector(MWTree<D, T> &tree) const {
    auto *nodeVec = new MWNodeVector<D, T>;
    if (tree.isPeriodic()) touchParentNodes(tree);
    tree_utils::make_node_table(tree, *nodeVec);
    return nodeVec;
}

// Explicit instantiations
template class ConvolutionCalculator<1, double>;
template class ConvolutionCalculator<2, double>;
template class ConvolutionCalculator<3, double>;

template class ConvolutionCalculator<1, ComplexDouble>;
template class ConvolutionCalculator<2, ComplexDouble>;
template class ConvolutionCalculator<3, ComplexDouble>;

} // namespace mrcpp
