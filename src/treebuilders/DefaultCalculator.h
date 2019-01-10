#pragma once

#include "TreeCalculator.h"

namespace mrcpp {

template<int D>
class DefaultCalculator final : public TreeCalculator<D> {
public:
    // Reimplementation without OpenMP, the default is faster this way
    void calcNodeVector(MWNodeVector<D> &nodeVec) {
        int nNodes = nodeVec.size();
        for (int n = 0; n < nNodes; n++) {
            calcNode(*nodeVec[n]);
        }
    }
private:
    void calcNode(MWNode<D> &node) {
        node.clearHasCoefs();
        node.clearNorms();
    }
};

}
