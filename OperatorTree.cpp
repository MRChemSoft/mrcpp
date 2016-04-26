#include "OperatorTree.h"
#include "OperatorNode.h"
#include "BandWidth.h"
#include "LebesgueIterator.h"

using namespace Eigen;
using namespace std;

OperatorTree::OperatorTree(const MultiResolutionAnalysis<2> &mra)
        : MWTree<2>(mra),
          bandWidth(0),
          nodePtrStore(0),
          nodePtrAccess(0) {
    for (int rIdx = 0; rIdx < this->rootBox.size(); rIdx++) {
        const NodeIndex<2> &nIdx = this->rootBox.getNodeIndex(rIdx);
        MWNode<2> *root = new OperatorNode(*this, nIdx);
        this->rootBox.setNode(rIdx, &root);
    }
    this->resetEndNodeTable();
}

OperatorTree::~OperatorTree() {
    clearOperNodeCache();
    clearBandWidth();
    MWNode<2> **roots = this->rootBox.getNodes();
    for (int i = 0; i < this->rootBox.size(); i++) {
        OperatorNode *node = static_cast<OperatorNode *>(roots[i]);
        if (node != 0) delete node;
        roots[i] = 0;
    }
}

void OperatorTree::clearBandWidth() {
    if (this->bandWidth != 0) delete this->bandWidth;
    this->bandWidth = 0;
}

void OperatorTree::calcBandWidth(double prec) {
    if (this->bandWidth != 0) MSG_ERROR("Band width not properly cleared");
    this->bandWidth = new BandWidth(getDepth());
    LebesgueIterator<2> it(this);
    double thrs = -1.0;
    while (it.next()) {
        MWNode<2> &node = it.getNode();
        const int *l = node.getNodeIndex().getTranslation();
        if (l[1] == 0) {
            int depth = node.getDepth();
            if (prec > 0.0) {
                thrs = max(MachinePrec, prec/(8.0 * (1 << depth)));
            }
            for (int k = 0; k < 4; k++) {
                if (node.getComponentNorm(k) > thrs) {
                    int bw = this->bandWidth->getWidth(depth, k);
                    if (l[0] > bw) {
                        this->bandWidth->setWidth(depth, k, l[0]);
                    }
                }
            }
        }
    }
}

void OperatorTree::getMaxTranslations(VectorXi &maxTransl) {
    int nScales = this->nodesAtDepth.size();
    maxTransl = VectorXi::Zero(nScales);
    LebesgueIterator<2> it(this);
    while(it.next()) {
        int n = it.getNode().getDepth();
        const int *l = it.getNode().getTranslation();
        maxTransl[n] = max(maxTransl[n], abs(l[0]));
        maxTransl[n] = max(maxTransl[n], abs(l[1]));
    }
}

/** Make 1D lists, adressable from [-l, l] scale by scale, of operator node
  * pointers for fast operator retrieval. This method is not thread safe,
  * since it projects missing operator nodes on the fly. Hence, it must NEVER
  * be called within a parallel region, or all hell will break loose. This is
  * not really a problem, but you have been warned.
  */
void OperatorTree::setupOperNodeCache() {
    int nScales = this->nodesAtDepth.size();
    int rootScale = this->getRootScale();
    this->nodePtrStore = new OperatorNode **[nScales];
    this->nodePtrAccess = new OperatorNode **[nScales];
    VectorXi max_transl;
    getMaxTranslations(max_transl);
    for (int n = 0; n < nScales; n++) {
        int scale = rootScale + n;
        int n_transl = max_transl[n];
        int n_nodes = 2*n_transl + 1;

        OperatorNode **nodes = new OperatorNode *[n_nodes];
        int j = 0;
        for (int i = n_transl; i >= 0; i--) {
            int l[2] = {0, i};
            NodeIndex<2> idx(scale, l);
            // Generated OperatorNodes are still OperatorNodes
            if (OperatorNode *oNode =
                dynamic_cast<OperatorNode *>(&MWTree<2>::getNode(idx))) {
                nodes[j] = oNode;
                j++;
            } else {
                NOT_REACHED_ABORT;
            }
        }
        for (int i = 1; i <= n_transl; i++) {
            int l[2] = {i, 0};
            NodeIndex<2> idx(scale, l);
            if (OperatorNode *oNode =
                dynamic_cast<OperatorNode *>(&MWTree<2>::getNode(idx))) {
                nodes[j] = oNode;
                j++;
            } else {
                NOT_REACHED_ABORT;
            }
        }

        this->nodePtrStore[n] = nodes;
        this->nodePtrAccess[n] = &nodes[n_transl];
    }
    this->resetEndNodeTable();
}

void OperatorTree::clearOperNodeCache() {
    if (this->nodePtrStore != 0 ) {
        for (int i = 0; i < getDepth(); i++) {
            delete[] this->nodePtrStore[i];
        }
        delete[] this->nodePtrStore;
        delete[] this->nodePtrAccess;
    }
}

