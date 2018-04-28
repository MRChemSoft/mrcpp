#include "treebuilders/grid.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/AnalyticAdaptor.h"
#include "treebuilders/CopyAdaptor.h"
#include "treebuilders/WaveletAdaptor.h"
#include "treebuilders/DefaultCalculator.h"
#include "mwutils/Printer.h"

namespace mrcpp {

template<int D>
void build_grid(FunctionTree<D> &out,
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
void copy_grid(FunctionTree<D> &out,
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
void copy_grid(FunctionTree<D> &out,
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
int clear_grid(double prec, FunctionTree<D> &out) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    DefaultCalculator<D> calculator;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    return builder.clear(out, calculator, adaptor);
}

template void build_grid(FunctionTree<1> &out, const RepresentableFunction<1> &inp, int maxIter);
template void build_grid(FunctionTree<2> &out, const RepresentableFunction<2> &inp, int maxIter);
template void build_grid(FunctionTree<3> &out, const RepresentableFunction<3> &inp, int maxIter);
template void copy_grid(FunctionTree<1> &out, FunctionTree<1> &inp, int maxIter);
template void copy_grid(FunctionTree<2> &out, FunctionTree<2> &inp, int maxIter);
template void copy_grid(FunctionTree<3> &out, FunctionTree<3> &inp, int maxIter);
template void copy_grid(FunctionTree<1> &out, FunctionTreeVector<1> &inp, int maxIter);
template void copy_grid(FunctionTree<2> &out, FunctionTreeVector<2> &inp, int maxIter);
template void copy_grid(FunctionTree<3> &out, FunctionTreeVector<3> &inp, int maxIter);
template int clear_grid(double prec, FunctionTree<1> &out);
template int clear_grid(double prec, FunctionTree<2> &out);
template int clear_grid(double prec, FunctionTree<3> &out);

} //namespace mrcpp
