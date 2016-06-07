#ifndef IDENTITYOPERATOR_H
#define IDENTITYOPERATOR_H

#include "ConvolutionOperator.h"
#include "IdentityKernel.h"

template<int D>
class IdentityOperator : public ConvolutionOperator<D> {
public:
    IdentityOperator(const MultiResolutionAnalysis<D> &mra,
                     double apply = -1.0, double build = -1.0)
            : ConvolutionOperator<D>(mra, apply, build, -1) {
        double epsilon = this->build_prec/10.0;
        IdentityKernel identity_kernel(epsilon);
        this->initializeOperator(identity_kernel);
    }
    virtual ~IdentityOperator() { }
};

#endif // IDENTITYOPERATOR_H
