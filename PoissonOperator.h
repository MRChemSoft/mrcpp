#ifndef POISSONOPERATOR_H
#define POISSONOPERATOR_H

#include "ConvolutionOperator.h"
#include "PoissonKernel.h"

class PoissonOperator : public ConvolutionOperator<3> {
public:
    PoissonOperator(const MultiResolutionAnalysis<3> &mra,
                    double apply = -1.0,
                    double build = -1.0,
                    int iter = -1)
            : ConvolutionOperator<3>(mra, apply, build, iter) {
        double epsilon = this->build_prec/100.0;
        double r_min = this->build_prec;
        double r_max = sqrt(3.0);
        PoissonKernel poisson_kernel(epsilon, r_min, r_max);
        initializeOperator(poisson_kernel);
    }
    virtual ~PoissonOperator() { }
};

#endif // POISSONOPERATOR_H
