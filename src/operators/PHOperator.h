#pragma once

#include "DerivativeOperator.h"

namespace mrcpp {

template<int D>
class PHOperator final : public DerivativeOperator<D> {
public:
    PHOperator(const MultiResolutionAnalysis<D> &mra, int order);
    PHOperator(const PHOperator &oper) = delete;
    PHOperator &operator=(const PHOperator &oper) = delete;

protected:
    void initializeOperator(int order);
};

}
