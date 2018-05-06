#pragma once

#include <functional>

#include "RepresentableFunction.h"

namespace mrcpp {

template<int D>
class AnalyticFunction : public RepresentableFunction<D> {
public:
    AnalyticFunction() { }
    AnalyticFunction(std::function<double (const double *r)> f,
                     const double *a = 0,
                     const double *b = 0)
        : RepresentableFunction<D>(a, b),
          func(f) { }
    virtual ~AnalyticFunction() { }

    void set(std::function<double (const double *r)> f) { this->func = f; }

    virtual double evalf(const double *r) const {
        double val = 0.0;
        if (not this->outOfBounds(r)) val = this->func(r);
        return val;
    }
protected:
    std::function<double (const double *r)> func;
};

}
