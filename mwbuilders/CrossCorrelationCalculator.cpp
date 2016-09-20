#include "CrossCorrelationCalculator.h"
#include "FunctionTree.h"
#include "MWNode.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

void CrossCorrelationCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();
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
    node.calcNorms();
}

template<int T>
void CrossCorrelationCalculator::applyCcc(MWNode<2> &node,
                                          CrossCorrelationCache<T> &ccc) {
    const MatrixXd &lMat = ccc.getLMatrix(node.getOrder());
    const MatrixXd &rMat = ccc.getRMatrix(node.getOrder());

    int scale = node.getScale() + 1;
    int t_dim = node.getMWTree().getTDim();
    int kp1_d = node.getKp1_d();

    VectorXd vec_o = VectorXd::Zero(t_dim*kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();
    for (int i = 0; i < t_dim; i++) {
        NodeIndex<2> cIdx(idx, i);
        const int *l = cIdx.getTranslation();
        int l_a = l[1] - l[0] - 1;
        int l_b = l[1] - l[0];

        NodeIndex<1> idx_a(scale, &l_a);
        NodeIndex<1> idx_b(scale, &l_b);

        const MWNode<1> &node_a = this->kernel->getNode(idx_a);
        const MWNode<1> &node_b = this->kernel->getNode(idx_b);

        VectorXd vec_a;
        VectorXd vec_b;
        node_a.getCoefs(vec_a);
        node_b.getCoefs(vec_b);

        const VectorXd &seg_a = vec_a.segment(0, node_a.getKp1_d());
        const VectorXd &seg_b = vec_b.segment(0, node_b.getKp1_d());
        vec_o.segment(i*kp1_d, kp1_d) = (lMat*seg_a + rMat*seg_b);
    }
    double *coefs = node.getCoefs_d();
    double two_n = pow(2.0, -scale/2.0);
    for (int i = 0; i < t_dim*kp1_d; i++) {
        coefs[i] = two_n*vec_o(i);
    }
}
