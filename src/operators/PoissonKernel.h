#pragma once

#include "GreensKernel.h"

namespace mrcpp {

class PoissonKernel final : public GreensKernel {
public:
    PoissonKernel(double eps, double r_min, double r_max)
            : GreensKernel(eps, r_min, r_max) {
        initializeKernel();
    }

protected:
    void initializeKernel();
    std::ostream& print(std::ostream &o) const;
};

}
