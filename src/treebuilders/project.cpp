#include "treebuilders/project.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/WaveletAdaptor.h"
#include "treebuilders/ProjectionCalculator.h"
#include "functions/AnalyticFunction.h"
#include "trees/MultiResolutionAnalysis.h"
#include "trees/FunctionTree.h"
#include "mwutils/Printer.h"
#include "mwutils/Timer.h"

namespace mrcpp {

template<int D>
void project(double prec,
             FunctionTree<D> &out,
             std::function<double (const double *r)> func,
             int maxIter) {
    AnalyticFunction<D> inp(func);
    mrcpp::project(prec, out, inp, maxIter);
}

template<int D>
void project(double prec,
             FunctionTree<D> &out,
             RepresentableFunction<D> &inp,
             int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    ProjectionCalculator<D> calculator(inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Printer::printTime(10, "Time transform", trans_t);
    Printer::printSeparator(10, ' ');
}

template void project(double prec, FunctionTree<1> &out, RepresentableFunction<1> &inp, int maxIter);
template void project(double prec, FunctionTree<2> &out, RepresentableFunction<2> &inp, int maxIter);
template void project(double prec, FunctionTree<3> &out, RepresentableFunction<3> &inp, int maxIter);
template void project(double prec, FunctionTree<1> &out, std::function<double (const double *r)> func, int maxIter);
template void project(double prec, FunctionTree<2> &out, std::function<double (const double *r)> func, int maxIter);
template void project(double prec, FunctionTree<3> &out, std::function<double (const double *r)> func, int maxIter);

} //namespace mrcpp
