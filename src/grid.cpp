#include "grid.h"
#include "TreeBuilder.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"
#include "WaveletAdaptor.h"
#include "DefaultCalculator.h"
#include "Printer.h"

using namespace mrcpp;

template<int D>
void mrcpp::build_grid(FunctionTree<D> &out,
                       const RepresentableFunction<D> &inp,
                       int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    AnalyticAdaptor<D> adaptor(inp, maxScale);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    Printer::printSeparator(10, ' ');
}

template<int D>
void mrcpp::copy_grid(FunctionTree<D> &out,
                      FunctionTree<D> &inp,
                      int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, maxScale, 0);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    Printer::printSeparator(10, ' ');
}

template<int D>
void mrcpp::copy_grid(FunctionTree<D> &out,
                      FunctionTreeVector<D> &inp,
                      int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, maxScale, 0);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    Printer::printSeparator(10, ' ');
}

template<int D>
int mrcpp::clear_grid(double prec, FunctionTree<D> &out) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    DefaultCalculator<D> calculator;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    return builder.clear(out, calculator, adaptor);
}

template void mrcpp::build_grid(FunctionTree<1> &out, const RepresentableFunction<1> &inp, int maxIter);
template void mrcpp::build_grid(FunctionTree<2> &out, const RepresentableFunction<2> &inp, int maxIter);
template void mrcpp::build_grid(FunctionTree<3> &out, const RepresentableFunction<3> &inp, int maxIter);
template void mrcpp::copy_grid(FunctionTree<1> &out, FunctionTree<1> &inp, int maxIter);
template void mrcpp::copy_grid(FunctionTree<2> &out, FunctionTree<2> &inp, int maxIter);
template void mrcpp::copy_grid(FunctionTree<3> &out, FunctionTree<3> &inp, int maxIter);
template void mrcpp::copy_grid(FunctionTree<1> &out, FunctionTreeVector<1> &inp, int maxIter);
template void mrcpp::copy_grid(FunctionTree<2> &out, FunctionTreeVector<2> &inp, int maxIter);
template void mrcpp::copy_grid(FunctionTree<3> &out, FunctionTreeVector<3> &inp, int maxIter);
template int mrcpp::clear_grid(double prec, FunctionTree<1> &out);
template int mrcpp::clear_grid(double prec, FunctionTree<2> &out);
template int mrcpp::clear_grid(double prec, FunctionTree<3> &out);
