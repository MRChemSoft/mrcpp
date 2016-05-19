/*
 *
 *
 *  \date Jul 7, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef POISSONKERNEL_H_
#define POISSONKERNEL_H_

#include "GreensKernel.h"

class PoissonKernel: public GreensKernel {
public:
    PoissonKernel(double eps, double r_min, double r_max)
            : GreensKernel(eps, r_min, r_max) {
        initializeKernel();
    }
    virtual ~PoissonKernel() { }
protected:
    virtual void initializeKernel();
};

#endif /* POISSONKERNEL_H_ */
