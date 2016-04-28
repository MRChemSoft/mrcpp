#include "OperatorNode.h"
#include "MathUtils.h"

using namespace Eigen;

OperatorNode::OperatorNode(OperatorTree &t, const NodeIndex<2> &n)
        : MWNode<2>(t, n) {
    this->allocCoefs(this->getTDim());
    this->setIsEndNode();
}

OperatorNode::OperatorNode(OperatorNode &p, int c)
        : MWNode<2>(p, c) {
    this->allocCoefs(this->getTDim());
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
double OperatorNode::calcComponentNorm(int i, double thrs) const {
    if (thrs < 0.0) thrs = MachineZero;
    int kp1_d = this->getKp1_d();
    const VectorXd &coefs = getCoefs().segment(i*kp1_d, kp1_d);

    double norm = 0.0;
    double vectorNorm = coefs.norm();
    if (vectorNorm > thrs) {
        double normInf = MathUtils::matrixNormInfinity(coefs);
        double normOne = MathUtils::matrixNorm1(coefs);
        if (sqrt(normInf*normOne) > thrs) {
            norm = MathUtils::matrixNorm2(coefs);
            if (norm < thrs) {
                norm = 0.0;
            }
        }
    }
    return norm;
}
