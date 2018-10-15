#pragma once

#include "ConvolutionOperator.h"
#include "DerivativeKernel.h"

namespace mrcpp {

template<int D>
class DerivativeConvolution final : public ConvolutionOperator<D> {
public:
    DerivativeConvolution(int d, const MultiResolutionAnalysis<D> &mra, double apply = -1.0, double build = -1.0)
            : ConvolutionOperator<D>(mra, apply, build) {
        double epsilon = this->build_prec/10.0;
        DerivativeKernel derivative_kernel(epsilon);
        this->initializeOperator(derivative_kernel);
        this->setApplyDir(d);
    }
};

}
