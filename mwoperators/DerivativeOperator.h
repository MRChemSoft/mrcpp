#ifndef DERIVATIVEOPERATOR_H
#define DERIVATIVEOPERATOR_H

#include "MWOperator.h"

template<int D>
class DerivativeOperator : public MWOperator {
public:
    DerivativeOperator(const MultiResolutionAnalysis<D> &mra) 
            : MWOperator(mra.getOperatorMRA()) { }
    virtual ~DerivativeOperator() { }
};

#endif // DERIVATIVEOPERATOR_H
