#include "CrossCorrelationCalculator.h"
#include "FunctionTree.h"
#include "MWNode.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

void CrossCorrelationCalculator::calcNode(MWNode<2> &node) {
    int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
    case Interpol:
    {
        getCrossCorrelationCache(Interpol, ccc);
        applyCcc(node, ccc);
        break;
    }
    case Legendre:
    {
        getCrossCorrelationCache(Legendre, ccc);
        applyCcc(node, ccc);
        break;
    }
    default:
        MSG_ERROR("Invalid scaling type");
        break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    int depth = node.getDepth();
    double thrs = std::max(MachinePrec, (0.1*this->prec)/(8.0 * (1 << depth)));
    node.calcNorms(thrs);
}

template<int T>
void CrossCorrelationCalculator::applyCcc(MWNode<2> &node,
                                          CrossCorrelationCache<T> &ccc) {
    const MatrixXd &lMat = ccc.getLMatrix(node.getOrder());
    const MatrixXd &rMat = ccc.getRMatrix(node.getOrder());

    int scale = node.getScale() + 1;
    int tDim = node.getMWTree().getTDim();

    const NodeIndex<2> &idx = node.getNodeIndex();
    for (int i = 0; i < tDim; i++) {
        NodeIndex<2> cIdx(idx, i);
        const int *l = cIdx.getTranslation();
        int l1 = l[1] - l[0] - 1;
        int l2 = l[1] - l[0];

        NodeIndex<1> idx1(scale, &l1);
        NodeIndex<1> idx2(scale, &l2);

        const MWNode<1> &node1 = this->kernel->getNode(idx1);
        const MWNode<1> &node2 = this->kernel->getNode(idx2);

        const Eigen::VectorXd &a = node1.getCoefs().segment(0,node1.getKp1_d());
        const Eigen::VectorXd &b = node2.getCoefs().segment(0,node2.getKp1_d());

        int nCoefs = node.getKp1_d();
        node.getCoefs().segment(i*nCoefs, nCoefs) = (lMat*a + rMat*b);
    }
    double factor = pow(2.0, -scale/2.0);
    node.getCoefs() *= factor;
}
