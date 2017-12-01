#pragma once

#include "DerivativeOperator.h"

template<int D>
class PHOperator : public DerivativeOperator<D> {
public:
    PHOperator(const MultiResolutionAnalysis<D> &mra, int order);
    virtual ~PHOperator() { }

protected:
    void initializeOperator(int order);
};

