#include "MWOperator.h"
#include "BandWidth.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

template<int D>
void MWOperator<D>::clearOperator() {
    for (int i = 0; i < this->oper_exp.size(); i++) {
        if (this->oper_exp[i] != 0) delete this->oper_exp[i];
    }
    this->oper_exp.clear();
}

template<int D>
OperatorTree& MWOperator<D>::getComponent(int i) {
    if (this->oper_exp[i] == 0) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

template<int D>
const OperatorTree& MWOperator<D>::getComponent(int i) const {
    if (this->oper_exp[i] == 0) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

template<int D>
int MWOperator<D>::getMaxBandWidth(int depth) const {
    int maxWidth = -1;
    if (depth < 0 ) {
        maxWidth = this->bandMax.maxCoeff();
    } else if (depth < this->bandMax.size() ) {
        maxWidth = this->bandMax(depth);
    }
    return maxWidth;
}

template<int D>
void MWOperator<D>::clearBandWidths() {
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        this->oper_exp[i]->clearBandWidth();
    }
}

template<int D>
void MWOperator<D>::calcBandWidths(double prec) {
    int maxDepth = 0;
    // First compute BandWidths and find depth of the deepest component
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        OperatorTree &oTree = *this->oper_exp[i];
        oTree.calcBandWidth(prec);
        const BandWidth &bw = oTree.getBandWidth();
        int depth = bw.getDepth();
        if (depth > maxDepth) {
            maxDepth = depth;
        }
    }
    this->bandMax = VectorXi(maxDepth + 1);
    this->bandMax.setConstant(-1);
    // Find the largest effective bandwidth at each scale
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        const OperatorTree &oTree = *this->oper_exp[i];
        const BandWidth &bw = oTree.getBandWidth();
        for (int n = 0; n <= bw.getDepth(); n++) { // scale loop
            for (int j = 0; j < 4; j++) { //component loop
                int w = bw.getWidth(n, j);
                if (w > this->bandMax(n)) {
                    this->bandMax(n) = w;
                }
            }
        }
    }
    println(20, "  Maximum bandwidths:\n" << this->bandMax << std::endl);
}

template class MWOperator<1>;
template class MWOperator<2>;
template class MWOperator<3>;
