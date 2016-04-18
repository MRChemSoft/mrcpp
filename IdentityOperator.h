#ifndef IDENTITYOPERATOR_H
#define IDENTITYOPERATOR_H

#include "MWOperator.h"
#include "IdentityCalculator.h"
#include "FunctionTree.h"
#include "Timer.h"

template<int D>
class IdentityOperator : public MWOperator<D> {
public:
    IdentityOperator(const MultiResolutionAnalysis<D> &mra,
                double prec = -1.0, int iter = -1)
        : MWOperator<D>(mra, prec, iter) {}
    IdentityOperator(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
        : MWOperator<D>(mra, a, iter) {}
    virtual ~IdentityOperator() {
    }

protected:
    TreeCalculator<D> *initCalculator(FunctionTree<D> &inp) {
        return new IdentityCalculator<D>(inp);
    }
};

#endif // IDENTITYOPERATOR_H
