#ifndef DERIVATIVEOPERATOR_H
#define DERIVATIVEOPERATOR_H

#include "MWOperator.h"

template<int D>
class DerivativeOperator : public MWOperator<D> {
public:
    DerivativeOperator(const MultiResolutionAnalysis<D> &mra,
                       double a = 0.5, double b = 0.5);
    virtual ~DerivativeOperator();

    void grad(FunctionTreeVector<D> &out, FunctionTree<D> &inp);
    void div(FunctionTree<D> &out, FunctionTreeVector<D> &inp);

protected:
    const double A;
    const double B;

    void initializeOperator();
};

#endif // DERIVATIVEOPERATOR_H
