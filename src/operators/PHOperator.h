#pragma once

#include "DerivativeOperator.h"

namespace mrcpp {

template<int D>
class PHOperator final : public DerivativeOperator<D> {
public:
    PHOperator(const MultiResolutionAnalysis<D> &mra, int order);

protected:
    void initializeOperator(int order);
};

}
