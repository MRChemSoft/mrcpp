#ifndef HELMHOLTZOPERATOR_H
#define HELMHOLTZOPERATOR_H

#include "ConvolutionOperator.h"
#include "HelmholtzKernel.h"

class HelmholtzOperator : public ConvolutionOperator<3> {
public:
    HelmholtzOperator(const MultiResolutionAnalysis<3> &mra,
                      double m,
                      double apply = -1.0,
                      double build = -1.0,
                      int iter = -1)
            : ConvolutionOperator<3>(mra, apply, build, iter),
              mu(m) {
        int plevel = TelePrompter::getPrintLevel();
        TelePrompter::setPrintLevel(0);
        double epsilon = this->build_prec/10.0;
        double r_min = calcMinDistance(epsilon);
        double r_max = calcMaxDistance();
        HelmholtzKernel helmholtz_kernel(this->mu, epsilon, r_min, r_max);
        // Rescale for application in 3D
        helmholtz_kernel.rescale(3);
        initializeOperator(helmholtz_kernel);
        TelePrompter::setPrintLevel(plevel);
    }
    virtual ~HelmholtzOperator() { }

    double getMu() const { return this->mu; }
protected:
    const double mu;
};


#endif // HELMHOLTZOPERATOR_H
