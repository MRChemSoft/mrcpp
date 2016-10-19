#include "OperatorNode.h"
#include "MathUtils.h"

using namespace Eigen;
using namespace std;


/** Calculate one specific component norm of the OperatorNode.
  *
  * OperatorNorms are defined as matrix 2-norms that are expensive to calculate.
  * Thus we calculate some cheaper upper bounds for this norm for thresholding.
  * First a simple vector norm, then a product of the 1- and infinity-norm. */
double OperatorNode::calcComponentNorm(int i) const {
    int depth = getDepth();
    double prec = getOperTree().getNormPrecision();
    double thrs = max(MachinePrec, prec/(8.0 * (1 << depth)));

    VectorXd coef_vec;
    this->getCoefs(coef_vec);

    int kp1_d = this->getKp1_d();
    const VectorXd &comp_vec = coef_vec.segment(i*kp1_d, kp1_d);

    double norm = 0.0;
    double vecNorm = comp_vec.norm();
    if (vecNorm > thrs) {
        double infNorm = MathUtils::matrixNormInfinity(comp_vec);
        double oneNorm = MathUtils::matrixNorm1(comp_vec);
        if (sqrt(infNorm*oneNorm) > thrs) {
            double twoNorm = MathUtils::matrixNorm2(comp_vec);
            if (twoNorm > thrs) {
                norm = twoNorm;
            }
        }
    }
    return norm;
}
