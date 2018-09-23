#include "grid.h"
#include "add.h"
#include "TreeBuilder.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"
#include "SplitAdaptor.h"
#include "WaveletAdaptor.h"
#include "DefaultCalculator.h"
#include "utils/Printer.h"

namespace mrcpp {

/** @brief Build grid of based on info from analytic function
 *
 * @param[in,out] out Output tree to be built
 * @param[in] inp Input function
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The grid of the output function will be EXTENDED using the general algorithm:
 *  1) Loop through current leaf nodes of the output tree
 *  2) Refine node based on custom split check from the function
 *  3) Repeat until convergence or maxIter is reached
 *
 * This algorithm will start at whatever grid is present in the output tree when
 * the function is called. This algorithm requires that the functions
 *  isVisibleAtScale()
 *  isZeroOnInterval()
 * is implemented in the particular RepresentableFunction.
 *
 * A negative maxIter means no bound.
 *
 */
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

/** @brief Build grid based on another MW function representation
 *
 * @param[in,out] out Output tree to be built
 * @param[in] inp Input tree
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The grid of the output function will be EXTENDED with all existing nodes in
 * corresponding input function, using the general algorithm:
 *  1) Loop through current leaf nodes of the output tree
 *  2) Refine node if the corresponding node in the input has children
 *  3) Repeat until all input nodes are covered or maxIter is reached
 *
 * This algorithm will start at whatever grid is present in the output tree when
 * the function is called. This means that all nodes on the input tree will also
 * be in the final output tree, but NOT vice versa.
 *
 * A negative maxIter means no bound.
 *
 */
template<int D>
void build_grid(FunctionTree<D> &out,
                FunctionTree<D> &inp,
                int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, maxScale, 0);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    Printer::printSeparator(10, ' ');
}

/** @brief Build grid based on several MW function representation
 *
 * @param[in,out] out Output tree to be built
 * @param[in] inp Input tree vector
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The grid of the output function will be EXTENDED with all existing nodes in
 * all corresponding input functions, using the general algorithm:
 *  1) Loop through current leaf nodes of the output tree
 *  2) Refine node if the corresponding node in one of the inputs has children
 *  3) Repeat until all input nodes are covered or maxIter is reached
 *
 * This algorithm will start at whatever grid is present in the output tree when
 * the function is called. This means that the final output grid will contain
 * (at least) the union of the nodes of all input trees.
 *
 * A negative maxIter means no bound.
 *
 */
template<int D>
void build_grid(FunctionTree<D> &out,
                FunctionTreeVector<D> &inp,
                int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, maxScale, 0);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    Printer::printSeparator(10, ' ');
}

/** @brief Copy function from one tree onto the grid of another tree
 *
 * @param[in,out] out Output tree to be built
 * @param[in] inp Input tree
 *
 * This algorithm will start at whatever grid is present in the output tree when
 * the function is called:
 *  1) Loop through current leaf nodes of the output tree
 *  2) Copy MW coefs from the corresponding input node
 *
 */
template<int D>
void copy_func(FunctionTree<D> &out, FunctionTree<D> &inp) {
    FunctionTreeVector<D> tmp_vec;
    tmp_vec.push_back(std::make_tuple(1.0, &inp));
    add(-1.0, out, tmp_vec);
}

/** @brief Copy the grid of another MW function representation
 *
 * @param[in,out] out Output tree to be built
 * @param[in] inp Input tree
 *
 * The grid of the output function will be identical to the grid of the input
 * function, but without MW coefficients.
 *
 */
template<int D>
void copy_grid(FunctionTree<D> &out, FunctionTree<D> &inp) {
    out.clear();
    build_grid(out, inp);
}

/** @brief Clear the grid of a MW function representation
 *
 * @param[in,out] out Output function to be cleared
 *
 * This will clear all MW coefs in the existing nodes, thus leaving an empty
 * grid that can be reused by computing new MW coefs.
 *
 */
template<int D>
void clear_grid(FunctionTree<D> &out) {
    TreeBuilder<D> builder;
    DefaultCalculator<D> calculator;
    builder.clear(out, calculator);
}

/** @brief Refine the grid of a MW function representation
 *
 * @param[in,out] out Output tree to be refined
 * @param[in] scales Number of refinement levels
 *
 * This will split ALL nodes in the tree the given number of times, then it will
 * compute scaling coefs of the new nodes, thus leaving the function representation
 * unchanged, but on a larger grid.
 *
 */
template<int D>
int refine_grid(FunctionTree<D> &out, int scales) {
    int nSplit = 0;
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    SplitAdaptor<D> adaptor(maxScale, true); // Splits all nodes
    for (int n = 0; n < scales; n++) {
        nSplit += builder.split(out, adaptor, true); // Transfers coefs to children
    }
    return nSplit;
}

/** @brief Refine the grid of a MW function representation
 *
 * @param[in,out] out Output tree to be refined
 * @param[in] prec Precision for initial split check
 *
 * This will first perform a split check on the existing end nodes in the tree
 * based on the provided precision parameter, then it will compute scaling coefs
 * of the new nodes, thus leaving the function representation unchanged, but on a
 * larger grid.
 *
 */
template<int D>
int refine_grid(FunctionTree<D> &out, double prec) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    int nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

/** @brief Refine the grid of a MW function representation
 *
 * @param[in,out] out Output tree to be refined
 * @param[in] inp Input tree
 *
 * This will first perform a split check on the existing end nodes in the output
 * tree based on the structure of the input tree (same as build_grid), then it will
 * compute scaling coefs of the new nodes, thus leaving the function representation
 * unchanged, but on a larger grid.
 *
 */
template<int D> int refine_grid(FunctionTree<D> &out, FunctionTree<D> &inp) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, maxScale, 0);
    int nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

template void build_grid(FunctionTree<1> &out, const RepresentableFunction<1> &inp, int maxIter);
template void build_grid(FunctionTree<2> &out, const RepresentableFunction<2> &inp, int maxIter);
template void build_grid(FunctionTree<3> &out, const RepresentableFunction<3> &inp, int maxIter);
template void build_grid(FunctionTree<1> &out, FunctionTree<1> &inp, int maxIter);
template void build_grid(FunctionTree<2> &out, FunctionTree<2> &inp, int maxIter);
template void build_grid(FunctionTree<3> &out, FunctionTree<3> &inp, int maxIter);
template void build_grid(FunctionTree<1> &out, FunctionTreeVector<1> &inp, int maxIter);
template void build_grid(FunctionTree<2> &out, FunctionTreeVector<2> &inp, int maxIter);
template void build_grid(FunctionTree<3> &out, FunctionTreeVector<3> &inp, int maxIter);
template void copy_func(FunctionTree<1> &out, FunctionTree<1> &inp);
template void copy_func(FunctionTree<2> &out, FunctionTree<2> &inp);
template void copy_func(FunctionTree<3> &out, FunctionTree<3> &inp);
template void copy_grid(FunctionTree<1> &out, FunctionTree<1> &inp);
template void copy_grid(FunctionTree<2> &out, FunctionTree<2> &inp);
template void copy_grid(FunctionTree<3> &out, FunctionTree<3> &inp);
template void clear_grid(FunctionTree<1> &out);
template void clear_grid(FunctionTree<2> &out);
template void clear_grid(FunctionTree<3> &out);
template int refine_grid(FunctionTree<1> &out, int scales);
template int refine_grid(FunctionTree<2> &out, int scales);
template int refine_grid(FunctionTree<3> &out, int scales);
template int refine_grid(FunctionTree<1> &out, double prec);
template int refine_grid(FunctionTree<2> &out, double prec);
template int refine_grid(FunctionTree<3> &out, double prec);
template int refine_grid(FunctionTree<1> &out, FunctionTree<1> &inp);
template int refine_grid(FunctionTree<2> &out, FunctionTree<2> &inp);
template int refine_grid(FunctionTree<3> &out, FunctionTree<3> &inp);

} // namespace mrcpp
