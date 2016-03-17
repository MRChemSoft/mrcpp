#ifndef MWMULTIPLIER_H
#define MWMULTIPLIER_H

#include "TreeBuilder.h"
#include "MultiplicationCalculator.h"
#include "GridGenerator.h"

template<int D>
class MWMultiplier : public TreeBuilder<D> {
public:
    MWMultiplier(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        NOT_IMPLEMENTED_ABORT;
    }
    MWMultiplier(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        NOT_IMPLEMENTED_ABORT;
    }
    virtual ~MWMultiplier() {
        NOT_IMPLEMENTED_ABORT;
    }

    FunctionTree<D>* operator()(MultiplicationVector<D> &inp) {
        NOT_IMPLEMENTED_ABORT;
    }

    void operator()(FunctionTree<D> &out, MultiplicationVector<D> &inp) {
        NOT_IMPLEMENTED_ABORT;
    }
};


#endif // MWMULTIPLIER_H
