#include "apply.h"
#include "grid.h"
#include "add.h"
#include "TreeBuilder.h"
#include "CopyAdaptor.h"
#include "SplitAdaptor.h"
#include "WaveletAdaptor.h"
#include "DefaultCalculator.h"
#include "DerivativeCalculator.h"
#include "ConvolutionCalculator.h"
#include "operators/DerivativeOperator.h"
#include "operators/ConvolutionOperator.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @brief Application of MW integral convolution operator
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] oper Convolution operator to apply
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

/** @brief Application of MW derivative operator
 *
 * @param[in,out] out Output function to be built
 * @param[in] oper Derivative operator to apply
 * @param[in] inp Input function
 * @param[in] dir Direction of derivative
 *
 * The output function will be computed on a FIXED grid that is predetermined
 * by the type of derivative operator. For a strictly local operator (ABGV_00),
 * the grid is an exact copy of the input function. For operators that involve
 * also neighboring nodes (ABGV_55, PH) the base grid will be WIDENED by one node
 * in the direction of application (on each side).
 *
 * The input function should contain only empty root nodes.
 *
 */
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
    SplitAdaptor<D> apply_adaptor(maxScale, false); // Splits no nodes
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

/** @brief Calculation of gradient vector of a function
 *
 * @param[in] oper Derivative operator to apply
 * @param[in] inp Input function
 * @returns FunctionTreeVector containing the gradient
 *
 * The derivative operator is applied in each Cartesian direction to the
 * input function and appended to the output vector.
 *
 * The length of the output vector will be the template dimension D.
 *
 */
template<int D>
FunctionTreeVector<D> gradient(DerivativeOperator<D> &oper, FunctionTree<D> &inp) {
    FunctionTreeVector<D> out;
    for (int d = 0; d < D; d++) {
        FunctionTree<D> *grad_d = new FunctionTree<D>(inp.getMRA());
        apply(*grad_d, oper, inp, d);
        out.push_back(std::make_tuple(1.0, grad_d));
    }
    return out;
}

/** @brief Calculation of divergence of a function vector
 *
 * @param[in,out] out Output function
 * @param[in] oper Derivative operator to apply
 * @param[in] inp Input function vector
 *
 * The derivative operator is applied in each Cartesian direction to the
 * corresponding components of the input vector and added up to the final
 * output. The grid of the output is fixed as the union of the component
 * grids (including any derivative widening, see derivative apply()).
 *
 * The length of the input vector must be the same as the template dimension D.
 *
 */
template<int D>
void divergence(FunctionTree<D> &out, DerivativeOperator<D> &oper, FunctionTreeVector<D> &inp) {
    if (inp.size() != D) MSG_FATAL("Dimension mismatch");

    FunctionTreeVector<D> tmp_vec;
    for (int d = 0; d < D; d++) {
        double coef_d = get_coef(inp, d);
        FunctionTree<D> &func_d = get_func(inp, d);
        FunctionTree<D> *out_d = new FunctionTree<D>(func_d.getMRA());
        apply(*out_d, oper, func_d, d);
        tmp_vec.push_back(std::make_tuple(coef_d, out_d));
    }
    build_grid(out, tmp_vec);
    add(-1.0, out, tmp_vec, 0); // Addition on union grid
    clear(tmp_vec, true);
}

template void apply(double prec, FunctionTree<1> &out, ConvolutionOperator<1> &oper, FunctionTree<1> &inp, int maxIter);
template void apply(double prec, FunctionTree<2> &out, ConvolutionOperator<2> &oper, FunctionTree<2> &inp, int maxIter);
template void apply(double prec, FunctionTree<3> &out, ConvolutionOperator<3> &oper, FunctionTree<3> &inp, int maxIter);
template void apply(FunctionTree<1> &out, DerivativeOperator<1> &oper, FunctionTree<1> &inp, int dir);
template void apply(FunctionTree<2> &out, DerivativeOperator<2> &oper, FunctionTree<2> &inp, int dir);
template void apply(FunctionTree<3> &out, DerivativeOperator<3> &oper, FunctionTree<3> &inp, int dir);
template void divergence(FunctionTree<1> &out, DerivativeOperator<1> &oper, FunctionTreeVector<1> &inp);
template void divergence(FunctionTree<2> &out, DerivativeOperator<2> &oper, FunctionTreeVector<2> &inp);
template void divergence(FunctionTree<3> &out, DerivativeOperator<3> &oper, FunctionTreeVector<3> &inp);
template FunctionTreeVector<1> gradient(DerivativeOperator<1> &oper, FunctionTree<1> &inp);
template FunctionTreeVector<2> gradient(DerivativeOperator<2> &oper, FunctionTree<2> &inp);
template FunctionTreeVector<3> gradient(DerivativeOperator<3> &oper, FunctionTree<3> &inp);

} // namespace mrcpp
