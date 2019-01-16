#pragma once

#include "DerivativeOperator.h"

namespace mrcpp {

template<int D>
class BSOperator final : public DerivativeOperator<D> {
public:
    explicit BSOperator(const MultiResolutionAnalysis<D> &mra, int order);
    explicit BSOperator(const BSOperator &oper) = delete;
    BSOperator &operator=(const BSOperator &oper) = delete;

protected:
    void initializeOperator();
};

}
