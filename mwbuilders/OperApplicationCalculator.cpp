#include "OperApplicationCalculator.h"
#include "OperatorTreeVector.h"
#include "OperatorState.h"
#include "FunctionNode.h"
#include "OperatorNode.h"
#include "BandWidth.h"
#include "eigen_disable_warnings.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace std;
using namespace Eigen;

template<int D>
OperApplicationCalculator<D>::OperApplicationCalculator(int dir,
                                                        double p,
                                                        OperatorTreeVector &o,
                                                        FunctionTree<D> &f,
                                                        int depth)
        : applyDir(dir),
          maxDepth(depth),
          prec(p),
          oper(&o),
          fTree(&f) {
    if (this->maxDepth > MaxDepth) MSG_FATAL("Beyond MaxDepth");
    initBandSizes();
}

template<int D>
OperApplicationCalculator<D>::~OperApplicationCalculator() {
    this->operStat.flushNodeCounters();
    println(10, this->operStat);
    for (int i = 0; i < this->bandSizes.size(); i++) {
        delete this->bandSizes[i];
    }
}

/** Initialize the number of nodes formally within the bandwidth of an
 operator. The band size is used for thresholding. */
template<int D>
void OperApplicationCalculator<D>::initBandSizes() {
    for (int i = 0; i < this->oper->size(); i++) {
        const OperatorTree &oTree = this->oper->getComponent(i);
        const BandWidth &bw = oTree.getBandWidth();
        MatrixXi *bsize = new MatrixXi(this->maxDepth, this->nComp2 + 1);
        bsize->setZero();
        for (int j = 0; j < this->maxDepth; j++) {
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
    this->operStat.incrementGNodeCounters(gNode);

    // Get all nodes in f within the bandwith of O in g
    MWNodeVector *fBand = makeOperBand(gNode);

    MWTree<D> &gTree = gNode.getMWTree();
    double gThrs = gTree.getSquareNorm();
    if (gThrs > 0.0) {
        double nTerms = (double) this->oper->size();
        gThrs = this->prec*sqrt(gThrs/nTerms);
    }
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

    for (int d = 0; d < D; d++) {
        int oTransl = fTransl[d] - gTransl[d];

//  The following will check the actual band width in each direction.
//  Not needed if the thresholding at the end of this routine is active.
        int a = (os.gt & (1 << d)) >> d;
        int b = (os.ft & (1 << d)) >> d;
        int idx = (a << 1) + b;
        int w = oTree.getBandWidth().getWidth(depth, idx);
        if (abs(oTransl) > w) {
            return;
        }

        const OperatorNode &oNode = oTree.getNode(depth, oTransl);
        int oIdx = os.getOperIndex(d);
        double ocn = oNode.getComponentNorm(oIdx);
        oNorm *= ocn;
        if (this->applyDir < 0 or this->applyDir == d) {
            oData[d] = const_cast<double *>(oNode.getCoefs_d()) + oIdx*os.kp1_2;
        } else {
            if (oTransl == 0 and (oIdx == 0 or oIdx == 3)) {
                // This will activate the identity operator in direction i
                oData[d] = 0;
            } else {
                // This means that we are in a zero part of the identity operator
                return;
            }
        }
    }
    double upperBound = oNorm * os.fThreshold;
    if (upperBound > os.gThreshold) {
        this->operStat.incrementFNodeCounters(fNode, os.ft, os.gt);
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
        if (oData[i] != 0) {
            if (i == D - 1) { // Last dir: Add up into g
                mult = 1.0;
            }
            const double *f = aux[i];
            double *g = const_cast<double *>(aux[i + 1]);
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    os.kp1_dm1, os.kp1, os.kp1, 1.0, f,
                    os.kp1, oData[i], os.kp1, mult, g, os.kp1_dm1);
        } else {
            // Identity operator in direction i
            Map<MatrixXd> f(aux[i], os.kp1, os.kp1_dm1);
            Map<MatrixXd> g(aux[i + 1], os.kp1_dm1, os.kp1);
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
        Map<MatrixXd> f(aux[i], os.kp1, os.kp1_dm1);
        Map<MatrixXd> g(aux[i + 1], os.kp1_dm1, os.kp1);
        if (oData[i] != 0) {
            Eigen::Map<MatrixXd> op(oData[i], os.kp1, os.kp1);
            if (i == D - 1) { // Last dir: Add up into g
                g += f.transpose() * op;
            } else {
                g = f.transpose() * op;
            }
        } else {
            // Identity operator in direction i
            if (i == D - 1) { // Last dir: Add up into g
                g += f.transpose();
            } else {
                g = f.transpose();
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

