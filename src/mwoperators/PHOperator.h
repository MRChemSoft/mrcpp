#pragma once

#include "mwoperators/DerivativeOperator.h"

namespace mrcpp {

template<int D>
class PHOperator : public DerivativeOperator<D> {
public:
    PHOperator(const MultiResolutionAnalysis<D> &mra, int order);
    virtual ~PHOperator() { }

protected:
    void initializeOperator(int order);
};

}
