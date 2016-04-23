#ifndef IDENTITYKERNEL_H
#define IDENTITYKERNEL_H

#include "GaussExp.h"
#include "GaussFunc.h"

class IdentityKernel : public GaussExp<1> {
public:
    IdentityKernel(double prec) {
        double alpha = sqrt(1.0/prec);
        double coef = pow(alpha/pi, 1.0/2.0);
        double pos = 0.0;
        GaussFunc<1> gauss(alpha, coef, &pos);
        this->append(gauss);
    }
    virtual ~IdentityKernel() { }
};

#endif // IDENTITYKERNEL_H
