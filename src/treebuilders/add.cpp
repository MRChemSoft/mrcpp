#include <tuple>
#include <vector>

#include "add.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "AdditionCalculator.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @brief Addition of two MW function representations
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] a Numerical coefficient of function a
 * @param[in] inp_a Input function a
 * @param[in] b Numerical coefficient of function b
 * @param[in] inp_b Input function b
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed as the sum of the two input functions
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
void add(double prec, FunctionTree<D> &out,
         double a, FunctionTree<D> &inp_a,
         double b, FunctionTree<D> &inp_b,
         int maxIter) {
    FunctionTreeVector<D> tmp_vec;
    tmp_vec.push_back(std::make_tuple(a, &inp_a));
    tmp_vec.push_back(std::make_tuple(b, &inp_b));
    add(prec, out, tmp_vec, maxIter);
}

/** @brief Addition of several MW function representations
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp Vector of input function
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed as the sum of all the functions
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
void add(double prec,
         FunctionTree<D> &out,
         FunctionTreeVector<D> &inp,
         int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    AdditionCalculator<D> calculator(inp);

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

template void add(double prec, FunctionTree<1> &out, double a, FunctionTree<1> &tree_a, double b, FunctionTree<1> &tree_b, int maxIter);
template void add(double prec, FunctionTree<2> &out, double a, FunctionTree<2> &tree_a, double b, FunctionTree<2> &tree_b, int maxIter);
template void add(double prec, FunctionTree<3> &out, double a, FunctionTree<3> &tree_a, double b, FunctionTree<3> &tree_b, int maxIter);
template void add(double prec, FunctionTree<1> &out, FunctionTreeVector<1> &inp, int maxIter);
template void add(double prec, FunctionTree<2> &out, FunctionTreeVector<2> &inp, int maxIter);
template void add(double prec, FunctionTree<3> &out, FunctionTreeVector<3> &inp, int maxIter);

} // namespace mrcpp
