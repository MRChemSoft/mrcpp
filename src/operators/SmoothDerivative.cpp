#include "SmoothDerivative.h"

namespace mrcpp {

template<int D>
SmoothDerivative<D>::SmoothDerivative(const MultiResolutionAnalysis<D>& mra,
                                      double prec,
                                      double cut_off,
                                      int max_Jpower)
: ConvolutionOperator<D>(mra, mra.getRootScale(), -10) {
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    // Keep this alive for the lifetime of the operator (avoid dangling pointer).
    this->cross_correlation =
        new SchrodingerEvolution_CrossCorrelation(30,
                                                  mra.getOrder(),
                                                  mra.getScalingBasis().getScalingType());

    initialize(cut_off, max_Jpower);

    this->initOperExp(1);
    Printer::setPrintLevel(oldlevel);
}

template<int D>
void SmoothDerivative<D>::initialize(double cut_off, int max_Jpower) {
    // Mirror TimeEvolutionOperator style to avoid API/method mismatches
    const int N = 18;

    const double o_prec = this->build_prec;
    auto        o_mra   = this->getOperatorMRA();

    auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);

    // J^{(deriv)} power integrals per scale
    std::map<int, DerivativePowerIntegrals*> J;
    for (int n = 0; n <= N + 1; ++n) {
        J[n] = new DerivativePowerIntegrals(cut_off, n, max_Jpower);
    }

    // Calculator uses the derivative series (real values) + cross-corr matrices
    DerivativeCrossCorrelationCalculator calculator(J, this->cross_correlation);

    // Adaptive build, identical pattern to TimeEvolutionOperator to ensure
    // CornerOperatorTree methods exist (mwTransform/removeRoughScaleNoise/etc.)
    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale(), true);
    TreeBuilder<2>  builder;
    builder.build(*o_tree, calculator, adaptor, N);

    // Postprocess
    Timer trans_t;
    o_tree->mwTransform(BottomUp);
    o_tree->removeRoughScaleNoise();
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; ++n) delete J[n];
}

// Explicit instantiations, matching your usual pattern
template class SmoothDerivative<1>;
template class SmoothDerivative<2>;
template class SmoothDerivative<3>;

} // namespace mrcpp