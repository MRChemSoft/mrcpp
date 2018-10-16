#pragma once

#include "MWOperator.h"

namespace mrcpp {

template<int D>
class DerivativeOperator : public MWOperator {
public:
    DerivativeOperator(const MultiResolutionAnalysis<D> &mra) 
            : MWOperator(mra.getOperatorMRA()) { }
    DerivativeOperator(const DerivativeOperator &oper) = delete;
    DerivativeOperator &operator=(const DerivativeOperator &oper) = delete;
    virtual ~DerivativeOperator() = default;
};

}
