#include "MWFDOperator.h"
#include "TreeBuilder.h"
#include "MWFDCalculator.h"
#include "WaveletAdaptor.h"
#include "FunctionTree.h"
#include "Timer.h"
#include "constants.h"

template<int D>
MWFDOperator<D>::MWFDOperator(const MultiResolutionAnalysis<D> &mra,
                              int k, int n, double pr)
    : prec(pr),
      maxScale(mra.getMaxScale()),
      diff_order(k),
      approx_order(n) {
    if (this->approx_order < 0)
        this->approx_order = mra.getOrder();
    if (this->diff_order > this->approx_order)
        MSG_FATAL("Derivative order must be lower than approximation order");
    if (mra.getOrder() < this->approx_order)
        MSG_FATAL("Approximation order must be lower than polynomial order");
}

template<int D>
void MWFDOperator<D>::operator()(FunctionTree<D> &out,
                                 FunctionTree<D> &inp,
                                 int maxIter,
                                 int dir) const {
    Timer pre_t;
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
    MWFDCalculator<D> calculator(dir, this->diff_order, this->approx_order, inp);
    pre_t.stop();

    builder.build(out, calculator, adaptor, maxIter);

    Timer post_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    inp.deleteGenerated();
    post_t.stop();

    println(10, "Time pre operator   " << pre_t);
    println(10, "Time post operator  " << post_t);
    println(10, std::endl);
}

template class MWFDOperator<1>;
template class MWFDOperator<2>;
template class MWFDOperator<3>;
