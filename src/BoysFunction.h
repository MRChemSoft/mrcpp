#pragma once
#include "RepresentableFunction.h"
#include "MWProjector.h"
#include "MultiResolutionAnalysis.h"

class BoysFunction : public RepresentableFunction<1> {
public:
    BoysFunction(int n, double prec = 1.0e-10);
    virtual ~BoysFunction() { }

    double evalf(const double *r) const;

protected:
    const int order;
    MWProjector<1> Q;
    MultiResolutionAnalysis<1> MRA;
};


