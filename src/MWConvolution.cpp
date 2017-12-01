#include "MWConvolution.h"
#include "WaveletAdaptor.h"
#include "ConvolutionCalculator.h"
#include "ConvolutionOperator.h"
#include "TreeBuilder.h"
#include "FunctionTree.h"
#include "Printer.h"
#include "Timer.h"

template<int D>
void MWConvolution<D>::operator()(FunctionTree<D> &out,
                                  ConvolutionOperator<D> &oper,
                                  FunctionTree<D> &inp,
                                  int max_iter) const {
    Timer pre_t;
    oper.calcBandWidths(this->prec);
    WaveletAdaptor<D> adaptor(this->prec, this->maxScale, this->absPrec);
    ConvolutionCalculator<D> calculator(this->prec, oper, inp);
    pre_t.stop();

    TreeBuilder<D> builder;
    builder.build(out, calculator, adaptor, max_iter);

    Timer post_t;
    oper.clearBandWidths();
    out.mwTransform(TopDown, false); // add coarse scale contributions
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    inp.deleteGenerated();
    post_t.stop();

    println(10, "Time pre operator   " << pre_t);
    println(10, "Time post operator  " << post_t);
    println(10, std::endl);
}

template class MWConvolution<1>;
template class MWConvolution<2>;
template class MWConvolution<3>;
