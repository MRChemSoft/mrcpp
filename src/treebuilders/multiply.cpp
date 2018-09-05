#include "multiply.h"
#include "add.h"
#include "grid.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "PowerCalculator.h"
#include "SquareCalculator.h"
#include "MultiplicationCalculator.h"
#include "trees/HilbertIterator.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionNode.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @brief Multiplication of two MW function representations
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] c Numerical coefficient
 * @param[in] inp_a Input function a
 * @param[in] inp_b Input function b
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed as the product of the two input functions
 * (including the numerical coefficient), using the general algorithm:
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
void multiply(double prec,
              FunctionTree<D> &out,
              double c,
              FunctionTree<D> &inp_a,
              FunctionTree<D> &inp_b,
              int maxIter) {
    FunctionTreeVector<D> tmp_vec;
    tmp_vec.push_back(std::make_tuple(c, &inp_a));
    tmp_vec.push_back(std::make_tuple(1.0, &inp_b));
    multiply(prec, out, tmp_vec, maxIter);
}

/** @brief Multiplication of several MW function representations
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp Vector of input function
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed as the product of all the functions
 * in the input vector (including their numerical coefficients), using the
 * general algorithm:
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
void multiply(double prec,
              FunctionTree<D> &out,
              FunctionTreeVector<D> &inp,
              int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    MultiplicationCalculator<D> calculator(inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    for (int i = 0; i < inp.size(); i++) {
        FunctionTree<D> &tree = get_func(inp, i);
        tree.deleteGenerated();
    }
    clean_t.stop();

    Printer::printTime(10, "Time transform", trans_t);
    Printer::printTime(10, "Time cleaning", clean_t);
    Printer::printSeparator(10, ' ');
}

/** @brief Out-of-place square of MW function representations
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp Input function to square
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed as the square of the input function,
 * using the general algorithm:
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
void square(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    SquareCalculator<D> calculator(inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    inp.deleteGenerated();
    clean_t.stop();

    Printer::printTime(10, "Time transform", trans_t);
    Printer::printTime(10, "Time cleaning", clean_t);
    Printer::printSeparator(10, ' ');
}

/** @brief Out-of-place power of MW function representations
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp Input function to square
 * @param[in] pow Numerical power
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed as the input function raised to the
 * given power, using the general algorithm:
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
void power(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, double pow, int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    PowerCalculator<D> calculator(inp, pow);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    inp.deleteGenerated();
    clean_t.stop();

    Printer::printTime(10, "Time transform", trans_t);
    Printer::printTime(10, "Time cleaning", clean_t);
    Printer::printSeparator(10, ' ');
}

template<int D>
void map(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, RepresentableFunction<D> &func) {
    NOT_IMPLEMENTED_ABORT;
}

/** @brief Dot product of two MW function vectors
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp_a Input function vector
 * @param[in] inp_b Input function vector
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed as the dot product of the two input
 * vectors (including their numerical coefficients). The precision parameter
 * is used only in the multiplication part, the final addition will be on
 * the fixed union grid of the components.
 *
 * The length of the input vectors must be the same.
 *
 */
template<int D>
void dot(double prec, FunctionTree<D> &out, FunctionTreeVector<D> &inp_a, FunctionTreeVector<D> &inp_b, int maxIter) {
    if (inp_a.size() != inp_b.size()) MSG_FATAL("Input length mismatch");

    FunctionTreeVector<D> tmp_vec;
    for (int d = 0; d < inp_a.size(); d++) {
        double coef_a = get_coef(inp_a, d);
        double coef_b = get_coef(inp_b, d);
        FunctionTree<D> &tree_a = get_func(inp_a, d);
        FunctionTree<D> &tree_b = get_func(inp_b, d);
        if (out.getMRA() != tree_a.getMRA()) MSG_FATAL("Trees not compatible");
        if (out.getMRA() != tree_b.getMRA()) MSG_FATAL("Trees not compatible");
        FunctionTree<D> *out_d = new FunctionTree<D>(out.getMRA());
        build_grid(*out_d, out);
        multiply(prec, *out_d, 1.0, tree_a, tree_b, maxIter);
        tmp_vec.push_back(std::make_tuple(coef_a*coef_b, out_d));
    }
    build_grid(out, tmp_vec);
    add(-1.0, out, tmp_vec, 0);
    clear(tmp_vec, true);
}

/** @brief Dot product of two MW function representations
 *
 * @param[in] bra Bra side input function
 * @param[in] ket Ket side input function
 *
 * The dot product is computed with the trees in compressed form, e.i. scaling
 * coefs only on root nodes, wavelet coefs on all nodes. Since wavelet functions
 * are orthonormal through ALL scales and the root scaling functions are
 * orthonormal to all finer level wavelet functions, this becomes a rather
 * efficient procedure as you only need to compute the dot product where the
 * grids overlaps.
 *
 */
template<int D>
double dot(FunctionTree<D> &bra, FunctionTree<D> &ket) {
    if (bra.getMRA() != ket.getMRA()){
        MSG_FATAL("Trees not compatible");
    }
    MWNodeVector nodeTable;
    HilbertIterator<D> it(&bra);
    it.setReturnGenNodes(false);
    while(it.next()) {
        MWNode<D> &node = it.getNode();
        nodeTable.push_back(&node);
    }
    int nNodes = nodeTable.size();
    double result = 0.0;
    double locResult = 0.0;
//OMP is disabled in order to get EXACT results (to the very last digit), the
//order of summation makes the result different beyond the 14th digit or so.
//OMP does improve the performace, but its not worth it for the time being.
//#pragma omp parallel firstprivate(n_nodes, locResult)
//		shared(nodeTable,rhs,result)
//    {
//#pragma omp for schedule(guided)
    for (int n = 0; n < nNodes; n++) {
        const FunctionNode<D> &braNode = static_cast<const FunctionNode<D> &>(*nodeTable[n]);
        const MWNode<D> *mwNode = ket.findNode(braNode.getNodeIndex());
        if (mwNode == 0) continue;

        const FunctionNode<D> &ketNode = static_cast<const FunctionNode<D> &>(*mwNode);
        if (braNode.isRootNode()) {
            locResult += dotScaling(braNode, ketNode);
        }
        locResult += dotWavelet(braNode, ketNode);
    }
//#pragma omp critical
    result += locResult;
//    }
    return result;
}

template void multiply(double prec, FunctionTree<1> &out, double c, FunctionTree<1> &tree_a, FunctionTree<1> &tree_b, int maxIter);
template void multiply(double prec, FunctionTree<2> &out, double c, FunctionTree<2> &tree_a, FunctionTree<2> &tree_b, int maxIter);
template void multiply(double prec, FunctionTree<3> &out, double c, FunctionTree<3> &tree_a, FunctionTree<3> &tree_b, int maxIter);
template void multiply(double prec, FunctionTree<1> &out, FunctionTreeVector<1> &inp, int maxIter);
template void multiply(double prec, FunctionTree<2> &out, FunctionTreeVector<2> &inp, int maxIter);
template void multiply(double prec, FunctionTree<3> &out, FunctionTreeVector<3> &inp, int maxIter);
template void power(double prec, FunctionTree<1> &out, FunctionTree<1> &tree, double pow, int maxIter);
template void power(double prec, FunctionTree<2> &out, FunctionTree<2> &tree, double pow, int maxIter);
template void power(double prec, FunctionTree<3> &out, FunctionTree<3> &tree, double pow, int maxIter);
template void square(double prec, FunctionTree<1> &out, FunctionTree<1> &tree, int maxIter);
template void square(double prec, FunctionTree<2> &out, FunctionTree<2> &tree, int maxIter);
template void square(double prec, FunctionTree<3> &out, FunctionTree<3> &tree, int maxIter);
template void map(double prec, FunctionTree<1> &out, FunctionTree<1> &inp, RepresentableFunction<1> &func);
template void map(double prec, FunctionTree<2> &out, FunctionTree<2> &inp, RepresentableFunction<2> &func);
template void map(double prec, FunctionTree<3> &out, FunctionTree<3> &inp, RepresentableFunction<3> &func);
template void dot(double prec, FunctionTree<1> &out, FunctionTreeVector<1> &inp_a, FunctionTreeVector<1> &inp_b, int maxIter);
template void dot(double prec, FunctionTree<2> &out, FunctionTreeVector<2> &inp_a, FunctionTreeVector<2> &inp_b, int maxIter);
template void dot(double prec, FunctionTree<3> &out, FunctionTreeVector<3> &inp_a, FunctionTreeVector<3> &inp_b, int maxIter);
template double dot(FunctionTree<1> &bra, FunctionTree<1> &ket);
template double dot(FunctionTree<2> &bra, FunctionTree<2> &ket);
template double dot(FunctionTree<3> &bra, FunctionTree<3> &ket);

} // namespace mrcpp
