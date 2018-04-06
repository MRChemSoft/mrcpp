/*
 * \breif
 */

#pragma once

#include "GaussExp.h"
#include "Gaussian.h"

namespace mrcpp {

class GreensKernel : public GaussExp<1> {
public:
    GreensKernel(double eps, double r_min, double r_max)
        : GaussExp<1>(),
          epsilon(eps),
          rMin(r_min),
          rMax(r_max) {
    }

    virtual ~GreensKernel() { }

    double getEpsilon() const { return this->epsilon; }
    double getRMin() const { return this->rMin; }
    double getRMax() const { return this->rMax; }

    void rescale(int d) {
        for (int i = 0; i < this->size(); i++) {
            Gaussian<1> &gauss = this->getFunc(i);
            double coef = pow(gauss.getCoef(), 1.0/d);
            gauss.setCoef(coef);
        }
    }

    friend std::ostream& operator <<(std::ostream &o, const GreensKernel &kernel) { return kernel.print(o); }

protected:
    const double epsilon;
    const double rMin; /**< lower extreme */
    const double rMax; /**< upper extreme  (not used) */

    virtual void initializeKernel() = 0;

    virtual std::ostream& print(std::ostream &o) const = 0;
};

}
