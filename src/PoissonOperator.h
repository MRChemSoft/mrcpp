#pragma once

#include "ConvolutionOperator.h"

class PoissonOperator : public ConvolutionOperator<3> {
public:
    PoissonOperator(const MultiResolutionAnalysis<3> &mra, double pr = -1.0);
    virtual ~PoissonOperator() { }
};

