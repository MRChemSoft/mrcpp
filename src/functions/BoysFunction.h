#pragma once

#include "RepresentableFunction.h"
#include "trees/MultiResolutionAnalysis.h"

namespace mrcpp {

class BoysFunction final : public RepresentableFunction<1> {
public:
    BoysFunction(int n, double prec = 1.0e-10);

    double evalf(const Coord<1> &r) const;

private:
    const int order;
    const double prec;
    MultiResolutionAnalysis<1> MRA;
};

}
