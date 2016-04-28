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
#include "macros.h"
#include "TelePrompter.h"

class PoissonKernel: public GreensKernel {
public:
    template<int D>
    PoissonKernel(const MultiResolutionAnalysis<D> &mra): GreensKernel(mra) {
        compKernRepr();
    }
    virtual ~PoissonKernel() {}
protected:
    virtual void compKernRepr();
};

#endif /* POISSONKERNEL_H_ */
