/*
 *
 *
 *  \date Jul 7, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef HELMHOLTZKERNEL_H_
#define HELMHOLTZKERNEL_H_

#include "GreensKernel.h"

class HelmholtzKernel: public GreensKernel {
public:
    HelmholtzKernel(double m, double eps, double r_min, double r_max)
            : GreensKernel(eps, r_min, r_max),
              mu(m) {
        initializeKernel();
    }
    virtual ~HelmholtzKernel() { }

    double getMu() const { return this->mu; }
protected:
    double mu; /**< exponent */
    virtual void initializeKernel();
};

#endif /* HELMHOLTZKERNEL_H_ */
