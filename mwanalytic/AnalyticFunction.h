#ifndef ANALYTICFUNCTION_H
#define ANALYTICFUNCTION_H

#include "RepresentableFunction.h"

template<int D>
class AnalyticFunction : public RepresentableFunction<D> {
public:
    AnalyticFunction(std::function<double(const double *r)> f,
                     const double *a = 0,
                     const double *b = 0)
        : RepresentableFunction<D>(a, b),
          func(f) { }
    virtual ~AnalyticFunction() { }

    virtual double evalf(const double *r) const {
        if (this->outOfBounds(r)) {
            return 0.0;
        } else {
            return this->func(r);
        }
    }
protected:
    std::function<double (const double *r)> func;
};

#endif // ANALYTICFUNCTION_H
