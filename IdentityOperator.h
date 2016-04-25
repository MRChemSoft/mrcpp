#ifndef IDENTITYOPERATOR_H
#define IDENTITYOPERATOR_H

#include "ConvolutionOperator.h"
#include "IdentityKernel.h"

template<int D>
class IdentityOperator : public ConvolutionOperator<D> {
public:
    IdentityOperator(const MultiResolutionAnalysis<D> &mra,
                     double apply = -1.0,
                     double build = -1.0,
                     int iter = -1)
            : ConvolutionOperator<D>(mra, apply, build, iter) {
        IdentityKernel identity_kernel(this->build_prec/100.0);
        this->initializeOperator(identity_kernel);
    }
    virtual ~IdentityOperator() { }
};

#endif // IDENTITYOPERATOR_H
