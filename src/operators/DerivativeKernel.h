#pragma once

#include "GreensKernel.h"
#include "functions/GaussFunc.h"
#include "functions/GaussPoly.h"

namespace mrcpp {

class DerivativeKernel final : public GreensKernel {
public:
    DerivativeKernel(double eps)
            : GreensKernel(eps, -1.0, -1.0) {
        initializeKernel();
    }
    DerivativeKernel(const DerivativeKernel &kern) = delete;
    DerivativeKernel &operator=(const DerivativeKernel &kern) = delete;

protected:
    void initializeKernel() {
        double alpha = 1.0/this->epsilon;
        double coef = std::pow(alpha/mrcpp::pi, 1.0/2.0);
        GaussFunc<1> g(alpha, coef);
        GaussPoly<1> dg = g.differentiate(0);
        this->append(dg);
    }

    std::ostream& print(std::ostream &o) const {
        o << " DerivativeKernel: " << std::endl;
        o << " epsilon:  " << this->epsilon << std::endl;
        return o;
    }
};

}
