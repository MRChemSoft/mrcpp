#include "OperatorNode.h"
#include "MathUtils.h"

using namespace Eigen;
using namespace std;

OperatorNode::OperatorNode(OperatorTree &t, const NodeIndex<2> &n)
        : MWNode<2>(t, n) {
    this->allocCoefs(this->getTDim(), this->getKp1_d());
    this->setIsEndNode();
}

OperatorNode::OperatorNode(OperatorNode &p, int c)
        : MWNode<2>(p, c) {
    this->allocCoefs(this->getTDim(), this->getKp1_d());
    this->setIsEndNode();
}

OperatorNode::~OperatorNode() {
    this->freeCoefs();
}

void OperatorNode::createChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MWNode<2> *child = new OperatorNode(*this, cIdx);
    this->children[cIdx] = child;
}

void OperatorNode::genChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MWNode<2> *child = new OperatorNode(*this, cIdx);
    this->children[cIdx] = child;
}

void OperatorNode::genChildren() {
    MWNode<2>::genChildren();
    this->giveChildrenCoefs();
}

/** Calculate one specific component norm of the OperatorNode.
  *
  * OperatorNorms are defined as matrix 2-norms that are expensive to calculate.
  * Thus we calculate some cheaper upper bounds for this norm for thresholding.
  * First a simple vector norm, then a product of the 1- and infinity-norm. */
double OperatorNode::calcComponentNorm(int i) const {
    NEEDS_FIX("Optimize without Eigen");
    int depth = getDepth();
    double prec = getOperTree().getNormPrecision();
    double thrs = max(MachinePrec, prec/(8.0 * (1 << depth)));

    int kp1_d = this->getKp1_d();
    VectorXd coefs;
    getCoefs(coefs);
    VectorXd comp_coefs = coefs.segment(i*kp1_d, kp1_d);

    double norm = 0.0;
    double vecNorm = comp_coefs.norm();
    if (vecNorm > thrs) {
        double infNorm = MathUtils::matrixNormInfinity(comp_coefs);
        double oneNorm = MathUtils::matrixNorm1(comp_coefs);
        if (sqrt(infNorm*oneNorm) > thrs) {
            double twoNorm = MathUtils::matrixNorm2(comp_coefs);
            if (twoNorm > thrs) {
                norm = twoNorm;
            }
        }
    }
    return norm;
}
