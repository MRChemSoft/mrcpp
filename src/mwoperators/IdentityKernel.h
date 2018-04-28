#pragma once

#include "mwoperators/GreensKernel.h"
#include "functions/GaussFunc.h"

namespace mrcpp {

class IdentityKernel : public GreensKernel {
public:
    IdentityKernel(double eps)
            : GreensKernel(eps, -1.0, -1.0) {
        initializeKernel();
    }
    virtual ~IdentityKernel() { }
protected:
    virtual void initializeKernel() {
        double alpha = sqrt(1.0/this->epsilon);
        double coef = pow(alpha/pi, 1.0/2.0);
        GaussFunc<1> gFunc(alpha, coef);
        this->append(gFunc);
    }
    virtual std::ostream& print(std::ostream &o) const {
        o << "Kernel: " << std::endl;
        o << "epsilon:  " << this->epsilon << std::endl;
        o << "rMin:     " << this->rMin << std::endl;
        o << "rMax:     " << this->rMax << std::endl;
        return o;
    }
};

}
