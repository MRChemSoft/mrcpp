#pragma once

#include "RepresentableFunction.h"
#include "trees/MultiResolutionAnalysis.h"

namespace mrcpp {

class BoysFunction : public RepresentableFunction<1> {
public:
    BoysFunction(int n, double prec = 1.0e-10);
    virtual ~BoysFunction() { }

    double evalf(const double *r) const;

protected:
    const int order;
    const double prec;
    MultiResolutionAnalysis<1> MRA;
};

}
