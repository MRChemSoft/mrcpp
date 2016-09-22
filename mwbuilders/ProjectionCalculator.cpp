#include <Eigen/Core>

#include "ProjectionCalculator.h"
#include "QuadratureCache.h"
#include "MultiResolutionAnalysis.h"
#include "MWNode.h"
#include "MWTree.h"

using namespace std;
using namespace Eigen;

template<int D>
void ProjectionCalculator<D>::calcNode(MWNode<D> &node) {
    const ScalingBasis &sf = node.getMWTree().getMRA().getScalingBasis();
    if (sf.getScalingType() != Interpol) {
        NOT_IMPLEMENTED_ABORT;
    }
    int quadratureOrder = sf.getQuadratureOrder();
    getQuadratureCache(qc);
    const VectorXd &pts = qc.getRoots(quadratureOrder);
    const VectorXd &wgts = qc.getWeights(quadratureOrder);

    double tmp_coefs[node.getNCoefs()];

    int scale = node.getScale();
    int kp1_d = node.getKp1_d();

    double scaleFactor = 1.0 / pow(2.0, scale + 1.0);
    double sqrtScaleFactor = sqrt(scaleFactor);
    double point[D];

    static int tDim = 1 << D;
    for (int cIdx = 0; cIdx < tDim; cIdx++) {
        NodeIndex<D> nIdx(node.getNodeIndex(), cIdx);
        const int *l = nIdx.getTranslation();

        int indexCounter[D];
        for (int i = 0; i < D; i++) {
            indexCounter[i] = 0;
        }

        for (int i = 0; i < kp1_d; i++) {
            double coef = 1.0;
            for (int j = 0; j < D; j++) {
                point[j] = scaleFactor * (pts(indexCounter[j]) + l[j]);
                coef *= sqrt(wgts(indexCounter[j])) * sqrtScaleFactor;
            }

            tmp_coefs[i] = coef * this->func->evalf(point);

            indexCounter[0]++;
            for (int j = 0; j < D - 1; j++) {
                if (indexCounter[j] == quadratureOrder) {
                    indexCounter[j] = 0;
                    indexCounter[j + 1]++;
                }
            }
        }
        node.setCoefBlock(cIdx, kp1_d, tmp_coefs);
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

template class ProjectionCalculator<1>;
template class ProjectionCalculator<2>;
template class ProjectionCalculator<3>;
