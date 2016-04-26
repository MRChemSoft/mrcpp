#include "OperApplicationCalculator.h"
#include "OperatorTreeVector.h"
#include "OperatorState.h"
#include "FunctionNode.h"
#include "OperatorNode.h"
#include "BandWidth.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace Eigen;

template<int D>
OperApplicationCalculator<D>::OperApplicationCalculator(OperatorTreeVector &o,
                                                        FunctionTree<D> &f)
        : oper(&o),
          fTree(&f) {
    initBandSizes();
    initNodeCounters();
}

template<int D>
OperApplicationCalculator<D>::~OperApplicationCalculator() {
    for (int i=0; i < this->bandSizes.size(); i++) {
        delete this->bandSizes[i];
    }
    for (int i = 0; i < this->nThreads; i++) {
        delete this->compCount[i];
        delete this->gNodeCount[i];
        delete this->fNodeCount[i];
    }
    delete [] this->compCount;
    delete [] this->fNodeCount;
    delete [] this->gNodeCount;
    delete [] this->genCount;
}

template<int D>
void OperApplicationCalculator<D>::initNodeCounters() {
    this->totAppNodes = 0;
    this->totGenAppNodes = 0;
    this->totGCount.setZero();
    this->totFCount.setZero();
    this->nThreads = omp_get_max_threads();
    this->compCount = new Eigen::Matrix<int, 8, 8> *[this->nThreads];
    this->fNodeCount = new Eigen::Vector2i *[this->nThreads];
    this->gNodeCount = new Eigen::Vector2i *[this->nThreads];
    this->genCount = new int[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->compCount[i] = new Eigen::Matrix<int, 8, 8>;
        this->compCount[i]->setZero();
        this->fNodeCount[i] = new Eigen::Vector2i;
        this->fNodeCount[i]->setZero();
        this->gNodeCount[i] = new Eigen::Vector2i;
        this->gNodeCount[i]->setZero();
        this->genCount[i] = 0;
    }
}

/** Reset all counters for operator application */
template<int D>
void OperApplicationCalculator<D>::resetNodeCounters() {
    NOT_IMPLEMENTED_ABORT;
    /*
    for (int i = 0; i < this->nThreads; i++) {
        this->fNodeCount[i]->setZero();
        this->gNodeCount[i]->setZero();
        this->compCount[i]->setZero();
        this->genCount[0] = 0;
    }
    this->totAppNodes = 0;
    this->totGenAppNodes = 0;
    this->totFCount.setZero();
    this->totGCount.setZero();
    this->fTree->clearNodeWeights();
    */
}

/** Sum all node counters from all threads. */
template<int D>
void OperApplicationCalculator<D>::flushNodeCounters() {
    NOT_IMPLEMENTED_ABORT;
    /*
    for (int i = 1; i < this->nThreads; i++) {
        *this->fNodeCount[0] += *this->fNodeCount[i];
        *this->gNodeCount[0] += *this->gNodeCount[i];
        *this->compCount[0] += *this->compCount[i];
        this->genCount[0] += this->genCount[i];
        this->compCount[i]->setZero();
        this->genCount[i] = 0;
    }
    this->totAppNodes = this->compCount[0]->sum();
    this->totGenAppNodes = this->genCount[0];
    this->totFCount = *this->fNodeCount[0];
    this->totGCount = *this->gNodeCount[0];
    */
}

/** Increment f-node usage counter. Needed for load balancing. */
template<int D>
void OperApplicationCalculator<D>::incrementFNodeCounters(MWNode<D> &gNode) {
    NOT_IMPLEMENTED_ABORT;
    /*
    int thread = omp_get_thread_num();
    if (gNode.isForeign()) {
        (*this->fNodeCount[thread])(1) += 1;
    } else {
        (*this->fNodeCount[thread])(0) += 1;
    }
    */
}

/** Increment g-node usage counter. Needed for load balancing. */
template<int D>
void OperApplicationCalculator<D>::incrementGNodeCounters(MWNode<D> &gNode) {
    NOT_IMPLEMENTED_ABORT;
    /*
    int thread = omp_get_thread_num();
    (*this->gNodeCount[thread])(1) += 1;
    */
}

/** Increment operator application counter. */
template<int D>
void OperApplicationCalculator<D>::incrementNodeCounters(int ft, int gt, bool isGenNode) {
    NOT_IMPLEMENTED_ABORT;
    /*
    int thread = omp_get_thread_num();
    (*this->compCount[thread])(ft,gt) += 1;
    if (isGenNode) {
        this->genCount[thread]++;
    }
    */
}

/** Initialize the number of nodes formally within the bandwidth of an
 operator. The band size is used for thresholding. */
template<int D>
void OperApplicationCalculator<D>::initBandSizes() {
    for (int i = 0; i < this->oper->size(); i++) {
        const OperatorTree &oTree = this->oper->getComponent(i);
        const BandWidth &bw = oTree.getBandWidth();
        MatrixXi *bsize = new MatrixXi(MaxDepth, this->nComp2 + 1);
        bsize->setZero();
        for (int j = 0; j < MaxDepth; j++) {
            calcBandSizeFactor(*bsize, j, bw);
        }
        this->bandSizes.push_back(bsize);
    }
}

/** Calculate the number of nodes within the bandwidth
  * of an operator. Currently this routine ignores the fact that
  * there are edges on the world box, and thus over estimates
  * the number of nodes. This is different from the previous version. */
template<int D>
void OperApplicationCalculator<D>::calcBandSizeFactor(MatrixXi &bs,
                                                      int depth,
                                                      const BandWidth &bw) {
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
template<int D>
MWNodeVector* OperApplicationCalculator<D>::makeOperBand(const MWNode<D> &gNode) {
    MWNodeVector *band = new MWNodeVector();

    int depth = gNode.getDepth();
    int width = this->oper->getMaxBandWidth(depth);
    if (width >= 0) {
        const NodeBox<D> &fWorld = this->fTree->getRootBox();
        const NodeIndex<D> &cIdx = fWorld.getCornerIndex();
        const NodeIndex<D> &gIdx = gNode.getNodeIndex();

        int l_start[D];
        int l_end[D];
        int nbox[D];
        for (int i = 0; i < D; i++) {
            l_start[i] = gIdx.getTranslation(i) - width;
            l_end[i] = gIdx.getTranslation(i) + width;
            // We need to consider the world borders
            int nboxes = fWorld.size(i) * (1 << depth);
            int c_i = cIdx.getTranslation(i) * (1 << depth);
            if (l_start[i] < c_i) {
                l_start[i] = c_i;
            }
            if (l_end[i] > c_i + nboxes - 1) {
                l_end[i] = c_i + nboxes - 1;
            }
            nbox[i] = l_end[i] - l_start[i] + 1;
        }

        int scale = gNode.getScale();
        NodeIndex<D> idx(scale, l_start);

        fillOperBand(band, idx, nbox, D-1);
    }
    return band;
}

/** Recursively retrieve all reachable f-nodes within the bandwidth. */
template<int D>
void OperApplicationCalculator<D>::fillOperBand(MWNodeVector *band,
                                                NodeIndex<D> &idx,
                                                const int *nbox,
                                                int dim) {
    int *l = idx.getTranslation();
    int l_start = l[dim];
    for (int j = 0; j < nbox[dim]; j++) {
        // Recurse until dim == 0
        if (dim > 0) {
            fillOperBand(band, idx, nbox, dim - 1);
            l[dim]++;
            continue;
        }
        MWNode<D> &fNode = this->fTree->getNode(idx);
        band->push_back(&fNode);
        l[dim]++;
    }
    l[dim] = l_start;
}

template<int D>
int OperApplicationCalculator<D>::getBandSizeFactor(int i, int depth,
                                                    const OperatorState<D> &os) const {
    assert(i >= 0 and i < this->bandSizes.size());
    MatrixXi &bs = *this->bandSizes[i];
    assert(depth >= 0 and depth <= bs.rows());
    int k = os.gt * this->nComp + os.ft;
    return bs(depth, k);
}

template<int D>
void OperApplicationCalculator<D>::calcNode(MWNode<D> &node) {
    FunctionNode<D> &gNode = static_cast<FunctionNode<D> &>(node);
    gNode.zeroCoefs();

    int depth = gNode.getDepth();
    OperatorState<D> os(gNode);
//    incrementGNodeCounters(gNode);

    // Get all nodes in f within the bandwith of O in g
    MWNodeVector *fBand = makeOperBand(gNode);

    MWTree<D> &gTree = gNode.getMWTree();
    double gThrs = gTree.getSquareNorm();
    //if (gThrs > 0.0) {
    //    gThrs = calcThreshold(gThrs, gNode.getMWTree().getRelPrec(),
    //            oper->getNTerms());
    //}
    os.gThreshold = gThrs;

    for (int n = 0; n < fBand->size(); n++) {
        MWNode<D> &fNode = *(*fBand)[n];
        os.setFNode(fNode);
        for (int ft = 0; ft < this->nComp; ft++) {
            double fNorm = fNode.getComponentNorm(ft);
            if (fNorm < MachineZero) {
                continue;
            }
            os.setFComponent(ft);
            for (int gt = 0; gt < this->nComp; gt++) {
                if (depth == 0 || gt != 0 || ft != 0) {
                    os.setGComponent(gt);
                    applyOperComp(os);
                }
            }
        }
    }
    gNode.calcNorms();
    delete fBand;
}

/** Apply each component (term) of the operator expansion to a node in f */
template<int D>
void OperApplicationCalculator<D>::applyOperComp(OperatorState<D> &os) {
    int depth = os.gNode->getDepth();
    double fNorm = os.fNode->getComponentNorm(os.ft);
    for (int i = 0; i < this->oper->size(); i++) {
        const OperatorTree &ot = this->oper->getComponent(i);
        const BandWidth &bw = ot.getBandWidth();
        if (os.getMaxDeltaL() > bw.getMaxWidth(depth)) {
            continue;
        }
        os.oTree = &ot;
        os.fThreshold = getBandSizeFactor(i, depth, os) * fNorm;
        applyOperator(os);
    }
}

/** Apply a single operator component (term) to a single f-node. Whether the
operator actualy is applied is determined by a screening threshold. */
template<int D>
void OperApplicationCalculator<D>::applyOperator(OperatorState<D> &os) {
    const OperatorTree &oTree = *os.oTree;
    MWNode<D> &gNode = *os.gNode;
    MWNode<D> &fNode = *os.fNode;

    const int *gTransl = gNode.getTranslation();
    const int *fTransl = fNode.getTranslation();
    int depth = gNode.getDepth();

    double oNorm = 1.0;
    double **oData = os.getOperData();

    for (int i = 0; i < D; i++) {
        int oTransl = fTransl[i] - gTransl[i];

//  The following will check the actual band width in each direction.
//  Not needed if the thresholding at the end of this routine is active.
        int a = (os.gt & (1 << i)) >> i;
        int b = (os.ft & (1 << i)) >> i;
        int idx = (a << 1) + b;
        int w = oTree.getBandWidth().getWidth(depth, idx);
        if (abs(oTransl) > w) {
            return;
        }

        const OperatorNode &oNode = oTree.getNode(depth, oTransl);
        int oIdx = os.getOperIndex(i);
        double ocn = oNode.getComponentNorm(oIdx);
        if (ocn == 0.0) { // Optimization. Not very relevant. Just a little.
            return;
        }
        oNorm *= ocn;
        oData[i] = const_cast<double *>(oNode.getCoefs().data()) + oIdx*os.kp1_2;
    }
    double upperBound = oNorm * os.fThreshold;
    if (upperBound > os.gThreshold) {
//        incrementNodeCounters(os.ft, os.gt, fNode.isGenNode());
//        incrementFNodeCounters(gNode);
        tensorApplyOperComp(os);
        gNode.setHasCoefs();
    }
}

/** Perorm the required linear algebra operations in order to apply an
operator component to a f-node in a n-dimensional tesor space. */
template<int D>
void OperApplicationCalculator<D>::tensorApplyOperComp(OperatorState<D> &os) {
    double **aux = os.getAuxData();
    double **oData = os.getOperData();
#ifdef HAVE_BLAS
    double mult = 0.0;
    for (int i = 0; i < D; i++) {
        if (oData[i] == 0) {
            Map<MatrixXd> f(aux[i], os.kp1, os.kp1_dm1);
            Map<MatrixXd> g(aux[i + 1], os.kp1_dm1, os.kp1);
            if (oData[i] == 0) {
                if (i == D - 1) { // Last dir: Add up into g
                    g += f.transpose();
                } else {
                    g = f.transpose();
                }
            }
        } else {
            if (i == D - 1) { // Last dir: Add up into g
                mult = 1.0;
            }
            const double *f = aux[i];
            double *g = const_cast<double *>(aux[i + 1]);
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    os.kp1_dm1, os.kp1, os.kp1, 1.0, f,
                    os.kp1, oData[i], os.kp1, mult, g, os.kp1_dm1);
        }
    }
#else
    for (int i = 0; i < D; i++) {
        Map<MatrixXd> f(aux[i], os.kp1, os.kp1_dm1);
        Map<MatrixXd> g(aux[i + 1], os.kp1_dm1, os.kp1);
        if (oData[i] == 0) {
            if (i == D - 1) { // Last dir: Add up into g
                g += f.transpose();
            } else {
                g = f.transpose();
            }
        } else {
            Eigen::Map<MatrixXd> op(oData[i], os.kp1, os.kp1);
            if (i == D - 1) { // Last dir: Add up into g
                g += f.transpose() * op;
            } else {
                g = f.transpose() * op;
            }
        }
    }
#endif
}

template<int D>
MWNodeVector* OperApplicationCalculator<D>::getInitialWorkVector(MWTree<D> &tree) const {
    MWNodeVector *nodeVec = new MWNodeVector;
    tree.makeNodeTable(*nodeVec);
    return nodeVec;
}

template class OperApplicationCalculator<1>;
template class OperApplicationCalculator<2>;
template class OperApplicationCalculator<3>;

