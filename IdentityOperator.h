#ifndef IDENTITYOPERATOR_H
#define IDENTITYOPERATOR_H

#include "MWOperator.h"

template<int D>
class IdentityOperator : public MWOperator<D> {
public:
    IdentityOperator(const MultiResolutionAnalysis<D> &mra,
                     double prec = -1.0, int iter = -1)
            : MWOperator<D>(mra, prec, iter) {
        initOperator();
    }
    virtual ~IdentityOperator() {
        this->oper.clear();
    }

protected:
    void initOperator() {
    }
};

#endif // IDENTITYOPERATOR_H
