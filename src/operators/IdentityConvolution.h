#pragma once

#include "ConvolutionOperator.h"
#include "IdentityKernel.h"

namespace mrcpp {

template<int D>
class IdentityConvolution final : public ConvolutionOperator<D> {
public:
    IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double pr = -1.0)
            : ConvolutionOperator<D>(mra, pr) {
        double epsilon = this->prec/10.0;
        IdentityKernel identity_kernel(epsilon);
        this->initializeOperator(identity_kernel);
    }
};

}
