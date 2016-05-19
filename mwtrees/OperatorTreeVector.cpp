#include "OperatorTreeVector.h"
#include "BandWidth.h"

using namespace std;
using namespace Eigen;

int OperatorTreeVector::getMaxBandWidth(int depth) const {
    int maxWidth = -1;
    if (depth < 0 ) {
        maxWidth = this->bandMax.maxCoeff();
    } else if (depth < this->bandMax.size() ) {
        maxWidth = this->bandMax(depth);
    }
    return maxWidth;
}

void OperatorTreeVector::clearBandWidths() {
    for (unsigned int i = 0; i < this->operComp.size(); i++) {
        this->operComp[i]->clearBandWidth();
    }
}

void OperatorTreeVector::calcBandWidths(double prec) {
    int maxDepth = 0;
    // First compute BandWidths and find depth of the deepest component
    for (unsigned int i = 0; i < this->operComp.size(); i++) {
        OperatorTree &oTree = *this->operComp[i];
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
    for (unsigned int i = 0; i < this->operComp.size(); i++) {
        const OperatorTree &oTree = *this->operComp[i];
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
