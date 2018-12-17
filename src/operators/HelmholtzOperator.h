#pragma once

#include "ConvolutionOperator.h"

namespace mrcpp {

class HelmholtzOperator final : public ConvolutionOperator<3> {
public:
    HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double m, double pr = -1.0);
    HelmholtzOperator(const HelmholtzOperator &oper) = delete;
    HelmholtzOperator &operator=(const HelmholtzOperator &oper) = delete;

    double getMu() const { return this->mu; }
protected:
    const double mu;
};

}
