#pragma once

#include "mwbuilders/TreeCalculator.h"
#include "mwcore/CrossCorrelationCache.h"

namespace mrcpp {

class CrossCorrelationCalculator : public TreeCalculator<2> {
public:
    CrossCorrelationCalculator(FunctionTree<1> &k) : kernel(&k) { }
    virtual ~CrossCorrelationCalculator() { }

protected:
    FunctionTree<1> *kernel;

    virtual void calcNode(MWNode<2> &node);

    template<int T>
    void applyCcc(MWNode<2> &node, CrossCorrelationCache<T> &ccc);
};

}
