#include "apply.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "ConvolutionCalculator.h"
#include "ConvolutionOperator.h"
#include "FunctionTree.h"
#include "Printer.h"
#include "Timer.h"

template<int D>
void mrcpp::apply(double prec,
                  FunctionTree<D> &out,
                  ConvolutionOperator<D> &oper,
                  FunctionTree<D> &inp,
                  int maxIter) {
    Timer pre_t;
    oper.calcBandWidths(prec);
    int maxScale = out.getMRA().getMaxScale();
    WaveletAdaptor<D> adaptor(prec, maxScale);
    ConvolutionCalculator<D> calculator(prec, oper, inp);
    pre_t.stop();

    TreeBuilder<D> builder;
    builder.build(out, calculator, adaptor, maxIter);

    Timer post_t;
    oper.clearBandWidths();
    out.mwTransform(TopDown, false); // add coarse scale contributions
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    inp.deleteGenerated();
    post_t.stop();

    Printer::printTime(10, "Time pre operator", pre_t);
    Printer::printTime(10, "Time post operator", post_t);
    Printer::printSeparator(10, ' ');
}

template void mrcpp::apply(double prec, FunctionTree<1> &out, ConvolutionOperator<1> &oper, FunctionTree<1> &inp, int maxIter);
template void mrcpp::apply(double prec, FunctionTree<2> &out, ConvolutionOperator<2> &oper, FunctionTree<2> &inp, int maxIter);
template void mrcpp::apply(double prec, FunctionTree<3> &out, ConvolutionOperator<3> &oper, FunctionTree<3> &inp, int maxIter);
