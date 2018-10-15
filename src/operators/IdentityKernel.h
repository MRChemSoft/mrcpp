#pragma once

#include "GreensKernel.h"
#include "functions/GaussFunc.h"

namespace mrcpp {

class IdentityKernel final : public GreensKernel {
public:
    IdentityKernel(double eps)
            : GreensKernel(eps, -1.0, -1.0) {
        initializeKernel();
    }
    IdentityKernel(const IdentityKernel &kern) = delete;
    IdentityKernel &operator=(const IdentityKernel &kern) = delete;

protected:
    void initializeKernel() {
        double alpha = std::sqrt(1.0/this->epsilon);
        double coef = std::pow(alpha/mrcpp::pi, 1.0/2.0);
        GaussFunc<1> gFunc(alpha, coef);
        this->append(gFunc);
    }
    std::ostream& print(std::ostream &o) const {
        o << "Kernel: " << std::endl;
        o << "epsilon:  " << this->epsilon << std::endl;
        o << "rMin:     " << this->rMin << std::endl;
        o << "rMax:     " << this->rMax << std::endl;
        return o;
    }
};

}
