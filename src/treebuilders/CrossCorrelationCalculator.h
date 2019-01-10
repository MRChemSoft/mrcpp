#pragma once

#include "TreeCalculator.h"
#include "core/CrossCorrelationCache.h"

namespace mrcpp {

class CrossCorrelationCalculator final : public TreeCalculator<2> {
public:
    CrossCorrelationCalculator(FunctionTree<1> &k) : kernel(&k) { }

private:
    FunctionTree<1> *kernel;

    void calcNode(MWNode<2> &node);

    template<int T>
    void applyCcc(MWNode<2> &node, CrossCorrelationCache<T> &ccc);
};

}
