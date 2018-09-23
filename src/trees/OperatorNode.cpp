#include "OperatorNode.h"
#include "SerialTree.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

void OperatorNode::dealloc() {
    int sIdx = this->serialIx;
    this->serialIx = -1;
    this->parentSerialIx = -1;
    this->childSerialIx = -1;
    this->tree->decrementNodeCount(this->getScale());
    this->tree->getSerialTree()->deallocNodes(sIdx);
}

/** Calculate one specific component norm of the OperatorNode.
  *
  * OperatorNorms are defined as matrix 2-norms that are expensive to calculate.
  * Thus we calculate some cheaper upper bounds for this norm for thresholding.
  * First a simple vector norm, then a product of the 1- and infinity-norm. */
double OperatorNode::calcComponentNorm(int i) const {
    int depth = getDepth();
    double prec = getOperTree().getNormPrecision();
    double thrs = std::max(MachinePrec, prec/(8.0 * (1 << depth)));

    VectorXd coef_vec;
    this->getCoefs(coef_vec);

    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    const VectorXd &comp_vec = coef_vec.segment(i*kp1_d, kp1_d);
    const MatrixXd comp_mat = MatrixXd::Map(comp_vec.data(), kp1, kp1);

    double norm = 0.0;
    double vecNorm = comp_vec.norm();
    if (vecNorm > thrs) {
        double infNorm = math_utils::matrix_norm_inf(comp_mat);
        double oneNorm = math_utils::matrix_norm_1(comp_mat);
        if (std::sqrt(infNorm*oneNorm) > thrs) {
            double twoNorm = math_utils::matrix_norm_2(comp_mat);
            if (twoNorm > thrs) norm = twoNorm;
        }
    }
    return norm;
}

} // namespace mrcpp
