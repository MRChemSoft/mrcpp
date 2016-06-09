#ifndef DERIVATIVEOPERATOR_H
#define DERIVATIVEOPERATOR_H

#include "ConvolutionOperator.h"
#include "DerivativeKernel.h"

template<int D>
class DerivativeOperator : public ConvolutionOperator<D> {
public:
    DerivativeOperator(const MultiResolutionAnalysis<D> &mra,
                       double apply = -1.0, double build = -1.0)
            : ConvolutionOperator<D>(mra, apply, build) {
        double epsilon = this->build_prec/10.0;
        DerivativeKernel derivative_kernel(epsilon);
        this->initializeOperator(derivative_kernel);
    }
    virtual ~DerivativeOperator() { }
};

#endif // DERIVATIVEOPERATOR_H
