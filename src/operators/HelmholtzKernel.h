#pragma once

#include "GreensKernel.h"

namespace mrcpp {

class HelmholtzKernel final : public GreensKernel {
public:
    HelmholtzKernel(double m, double eps, double r_min, double r_max)
            : GreensKernel(eps, r_min, r_max),
              mu(m) {
        initializeKernel();
    }
    HelmholtzKernel(const HelmholtzKernel &kern) = delete;
    HelmholtzKernel &operator=(const HelmholtzKernel &kern) = delete;

protected:
    const double mu; /**< exponent */
    void initializeKernel() override;
    std::ostream& print(std::ostream &o) const;
};

}
