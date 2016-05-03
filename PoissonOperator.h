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
        int plevel = TelePrompter::getPrintLevel();
        TelePrompter::setPrintLevel(0);
        double epsilon = this->build_prec/100.0;
        double r_min = calcMinDistance(epsilon);
        double r_max = calcMaxDistance();
        PoissonKernel poisson_kernel(epsilon, r_min, r_max);
        // Rescale for application in 3D
        poisson_kernel.rescale(3);
        initializeOperator(poisson_kernel);
        TelePrompter::setPrintLevel(plevel);
    }
    virtual ~PoissonOperator() { }
};

#endif // POISSONOPERATOR_H
