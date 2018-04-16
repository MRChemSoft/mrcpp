#include "mwbuilders/apply.h"
#include "mwbuilders/TreeBuilder.h"
#include "mwbuilders/CopyAdaptor.h"
#include "mwbuilders/WaveletAdaptor.h"
#include "mwbuilders/DefaultCalculator.h"
#include "mwbuilders/DerivativeCalculator.h"
#include "mwbuilders/ConvolutionCalculator.h"
#include "mwoperators/DerivativeOperator.h"
#include "mwoperators/ConvolutionOperator.h"
#include "mwtrees/FunctionTree.h"
#include "mwutils/Printer.h"
#include "mwutils/Timer.h"

namespace mrcpp {

template<int D>
void apply(double prec,
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

template<int D>
void apply(FunctionTree<D> &out,
           DerivativeOperator<D> &oper,
           FunctionTree<D> &inp,
           int dir) {
    TreeBuilder<D> builder;
    int maxScale = out.getMRA().getMaxScale();

    int bw[D]; // Operator bandwidth in [x,y,z]
    for (int d = 0; d < D; d++) bw[d] = 0;

    // Copy input tree plus bandwidth in operator direction
    Timer pre_t;
    oper.calcBandWidths(1.0); // Fixed 0 or 1 for derivatives
    bw[dir] = oper.getMaxBandWidth();
    CopyAdaptor<D> pre_adaptor(inp, maxScale, bw);
    DefaultCalculator<D> pre_calculator;
    builder.build(out, pre_calculator, pre_adaptor, -1);
    pre_t.stop();

    // Apply operator on fixed expanded grid
    TreeAdaptor<D> apply_adaptor(maxScale);
    DerivativeCalculator<D> apply_calculator(dir, oper, inp);
    builder.build(out, apply_calculator, apply_adaptor, 0);

    Timer post_t;
    oper.clearBandWidths();
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    inp.deleteGenerated();
    post_t.stop();

    Printer::printTime(10, "Time pre operator", pre_t);
    Printer::printTime(10, "Time post operator", post_t);
    Printer::printSeparator(10, ' ');
}

template void apply(double prec, FunctionTree<1> &out, ConvolutionOperator<1> &oper, FunctionTree<1> &inp, int maxIter);
template void apply(double prec, FunctionTree<2> &out, ConvolutionOperator<2> &oper, FunctionTree<2> &inp, int maxIter);
template void apply(double prec, FunctionTree<3> &out, ConvolutionOperator<3> &oper, FunctionTree<3> &inp, int maxIter);
template void apply(FunctionTree<1> &out, DerivativeOperator<1> &oper, FunctionTree<1> &inp, int dir);
template void apply(FunctionTree<2> &out, DerivativeOperator<2> &oper, FunctionTree<2> &inp, int dir);
template void apply(FunctionTree<3> &out, DerivativeOperator<3> &oper, FunctionTree<3> &inp, int dir);

} //namespace mrcpp
