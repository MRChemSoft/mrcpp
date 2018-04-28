#pragma once

#include "operators/GreensKernel.h"

namespace mrcpp {

class PoissonKernel: public GreensKernel {
public:
    PoissonKernel(double eps, double r_min, double r_max)
            : GreensKernel(eps, r_min, r_max) {
        initializeKernel();
    }
    virtual ~PoissonKernel() { }
protected:
    virtual void initializeKernel();
    virtual std::ostream& print(std::ostream &o) const;
};

}
