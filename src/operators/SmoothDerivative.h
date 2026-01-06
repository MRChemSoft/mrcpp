#pragma once

#include "ConvolutionOperator.h"
#include "MWOperator.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"
#include "functions/JpowerIntegrals.h"                       // DerivativePowerIntegrals
#include "treebuilders/TimeEvolution_CrossCorrelationCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "trees/CornerOperatorTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/**
 * SmoothDerivative: builds a smoothed first-derivative operator
 * (uses DerivativePowerIntegrals + DerivativeCrossCorrelationCalculator).
 *
 * Kept separate from TimeEvolutionOperator to avoid breaking your existing code.
 */
template<int D>
class SmoothDerivative : public ConvolutionOperator<D> {
public:
    SmoothDerivative(const MultiResolutionAnalysis<D>& mra,
                     double prec,
                     double cut_off,
                     int max_Jpower = 30);

    SmoothDerivative(const SmoothDerivative&)            = delete;
    SmoothDerivative& operator=(const SmoothDerivative&) = delete;
    ~SmoothDerivative() override                          = default;

    double getBuildPrec() const { return this->build_prec; }

protected:
    void initialize(double cut_off, int max_Jpower);
    void setBuildPrec(double prec) { this->build_prec = prec; }

    double build_prec{-1.0};
    SchrodingerEvolution_CrossCorrelation* cross_correlation{nullptr};
};

} // namespace mrcpp