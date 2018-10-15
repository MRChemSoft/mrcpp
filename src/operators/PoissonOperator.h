#pragma once

#include "ConvolutionOperator.h"

namespace mrcpp {

class PoissonOperator final : public ConvolutionOperator<3> {
public:
    PoissonOperator(const MultiResolutionAnalysis<3> &mra, double pr = -1.0);
};

}
