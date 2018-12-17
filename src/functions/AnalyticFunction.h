#pragma once

#include <functional>

#include "RepresentableFunction.h"

namespace mrcpp {

template<int D>
class AnalyticFunction : public RepresentableFunction<D> {
public:
    AnalyticFunction() = default;
    virtual ~AnalyticFunction() = default;

    AnalyticFunction(std::function<double (const Coord<D> &r)> f,
                     const double *a = nullptr,
                     const double *b = nullptr)
        : RepresentableFunction<D>(a, b), func(f) { }

    void set(std::function<double (const Coord<D> &r)> f) { this->func = f; }

    virtual double evalf(const Coord<D> &r) const {
        double val = 0.0;
        if (not this->outOfBounds(r)) val = this->func(r);
        return val;
    }
protected:
    std::function<double (const Coord<D> &r)> func;
};

}
