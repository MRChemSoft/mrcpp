#include "MWOperator.h"
#include "BandWidth.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

void MWOperator::clear(bool dealloc) {
    if (dealloc) {
        for (int i = 0; i < this->oper_exp.size(); i++) {
            if (this->oper_exp[i] != 0) delete this->oper_exp[i];
        }
    }
    this->oper_exp.clear();
}

OperatorTree& MWOperator::getComponent(int i) {
    if (this->oper_exp[i] == 0) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

const OperatorTree& MWOperator::getComponent(int i) const {
    if (this->oper_exp[i] == 0) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

int MWOperator::getMaxBandWidth(int depth) const {
    int maxWidth = -1;
    if (depth < 0 ) {
        maxWidth = this->band_max.maxCoeff();
    } else if (depth < this->band_max.size() ) {
        maxWidth = this->band_max(depth);
    }
    return maxWidth;
}

void MWOperator::clearBandWidths() {
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        this->oper_exp[i]->clearBandWidth();
    }
}

void MWOperator::calcBandWidths(double prec) {
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
    this->band_max = VectorXi(maxDepth + 1);
    this->band_max.setConstant(-1);
    // Find the largest effective bandwidth at each scale
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        const OperatorTree &oTree = *this->oper_exp[i];
        const BandWidth &bw = oTree.getBandWidth();
        for (int n = 0; n <= bw.getDepth(); n++) { // scale loop
            for (int j = 0; j < 4; j++) { //component loop
                int w = bw.getWidth(n, j);
                if (w > this->band_max(n)) {
                    this->band_max(n) = w;
                }
            }
        }
    }
    println(20, "  Maximum bandwidths:\n" << this->band_max << std::endl);
}
