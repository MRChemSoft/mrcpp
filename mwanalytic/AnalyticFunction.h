#ifndef ANALYTICFUNCTION_H
#define ANALYTICFUNCTION_H

#include "RepresentableFunction.h"

template<int D>
class AnalyticFunction : public RepresentableFunction<D> {
public:
    AnalyticFunction(double (*f)(const double *r), const double *a = 0, const double *b = 0)
        : RepresentableFunction<D>(a, b),
          func(f) { }
    virtual ~AnalyticFunction() { }

    virtual double evalf(const double *r) const {
        if (this->outOfBounds(r)) {
            return 0.0;
        } else {
            return (*this->func)(r);
        }
    }
protected:
    double (*func)(const double *r);
};

#endif // ANALYTICFUNCTION_H
