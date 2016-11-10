#ifndef IDENTITYOPERATOR_H
#define IDENTITYOPERATOR_H

#include "ConvolutionOperator.h"
#include "IdentityKernel.h"

template<int D>
class IdentityOperator : public ConvolutionOperator<D> {
public:
    IdentityOperator(const MultiResolutionAnalysis<D> &mra, double pr = -1.0)
            : ConvolutionOperator<D>(mra, pr) {
        double epsilon = this->prec/10.0;
        IdentityKernel identity_kernel(epsilon);
        this->initializeOperator(identity_kernel);
    }
    virtual ~IdentityOperator() { }
    bool applyCompressed() const { return false; }
};

#endif // IDENTITYOPERATOR_H
