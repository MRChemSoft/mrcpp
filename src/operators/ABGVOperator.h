#pragma once

#include "DerivativeOperator.h"

namespace mrcpp {

template<int D>
class ABGVOperator final : public DerivativeOperator<D> {
public:
    ABGVOperator(const MultiResolutionAnalysis<D> &mra, double a, double b);
    ABGVOperator(const ABGVOperator &oper) = delete;
    ABGVOperator &operator=(const ABGVOperator &oper) = delete;

protected:
    void initializeOperator(double a, double b);
};

} // namespace mrcpp
