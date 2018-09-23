#include "project.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "ProjectionCalculator.h"
#include "functions/AnalyticFunction.h"
#include "trees/MultiResolutionAnalysis.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @brief Projection of analytic function into MW representation
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp Input function
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed using the general algorithm:
 *  1) Compute MW coefs on current grid
 *  2) Refine grid where necessary based on prec
 *  3) Repeat until convergence or maxIter is reached
 *
 * This algorithm will start at whatever grid is present in the output tree when
 * the function is called (this grid should however be EMPTY, e.i. no coefs).
 *
 * A negative precision means NO refinement, as do maxIter = 0.
 * A negative maxIter means no bound.
 *
 */
template<int D>
void project(double prec,
             FunctionTree<D> &out,
             std::function<double (const double *r)> func,
             int maxIter) {
    AnalyticFunction<D> inp(func);
    mrcpp::project(prec, out, inp, maxIter);
}

/** @brief Projection of analytic function into MW representation
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp Input function
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed using the general algorithm:
 *  1) Compute MW coefs on current grid
 *  2) Refine grid where necessary based on prec
 *  3) Repeat until convergence or maxIter is reached
 *
 * This algorithm will start at whatever grid is present in the output tree when
 * the function is called (this grid should however be EMPTY, e.i. no coefs).
 *
 * A negative precision means NO refinement, as do maxIter = 0.
 * A negative maxIter means no bound.
 *
 */
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

} // namespace mrcpp
