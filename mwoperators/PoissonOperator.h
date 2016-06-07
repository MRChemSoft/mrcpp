#ifndef POISSONOPERATOR_H
#define POISSONOPERATOR_H

#include "ConvolutionOperator.h"
#include "PoissonKernel.h"

class PoissonOperator : public ConvolutionOperator<3> {
public:
    PoissonOperator(const MultiResolutionAnalysis<3> &mra,
                    double apply = -1.0,
                    double build = -1.0)
            : ConvolutionOperator<3>(mra, apply, build, -1) {
        int oldlevel = TelePrompter::setPrintLevel(0);
        double epsilon = this->build_prec/10.0;
        double r_min = calcMinDistance(epsilon);
        double r_max = calcMaxDistance();
        PoissonKernel poisson_kernel(epsilon, r_min, r_max);
        // Rescale for application in 3D
        poisson_kernel.rescale(3);
        initializeOperator(poisson_kernel);
        TelePrompter::setPrintLevel(oldlevel);
    }
    virtual ~PoissonOperator() { }
};

#endif // POISSONOPERATOR_H
