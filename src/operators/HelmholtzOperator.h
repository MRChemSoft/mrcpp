#pragma once

#include "operators/ConvolutionOperator.h"

namespace mrcpp {

class HelmholtzOperator : public ConvolutionOperator<3> {
public:
    HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double m, double pr = -1.0);
    virtual ~HelmholtzOperator() { }

    double getMu() const { return this->mu; }
protected:
    const double mu;
};

}
