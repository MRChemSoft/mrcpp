/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "grid.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"
#include "DefaultCalculator.h"
#include "SplitAdaptor.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "functions/GaussExp.h"
#include "functions/function_utils.h"
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
template <int D> void build_grid(FunctionTree<D> &out, const RepresentableFunction<D> &inp, int maxIter) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    AnalyticAdaptor<D> adaptor(inp, maxScale);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
}

template <int D> void build_grid(FunctionTree<D> &out, const Gaussian<D> &inp, int maxIter) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    DefaultCalculator<D> calculator;

    if (!out.getMRA().getWorldBox().isPeriodic()) {
        AnalyticAdaptor<D> adaptor(inp, maxScale);
        builder.build(out, calculator, adaptor, maxIter);
    } else {
        auto period = out.getMRA().getWorldBox().getScalingFactor();
        auto g_exp = function_utils::make_gaussian_periodic<D>(inp, period);
        for (auto i = 0; i < g_exp->size(); i++) {
            AnalyticAdaptor<D> adaptor(g_exp->getFunc(i), maxScale);
            builder.build(out, calculator, adaptor, maxIter);
        }
    }
    Printer::printSeparator(10, ' ');
}

template <int D> void build_grid(FunctionTree<D> &out, const GaussExp<D> &inp, int maxIter) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    DefaultCalculator<D> calculator;
    if (!out.getMRA().getWorldBox().isPeriodic()) {
        for (auto i = 0; i < inp.size(); i++) {
            AnalyticAdaptor<D> adaptor(inp.getFunc(i), maxScale);
            builder.build(out, calculator, adaptor, maxIter);
        }
    } else {
        auto period = out.getMRA().getWorldBox().getScalingFactor();
        auto *gauss = inp.getFunc(0).copy();
        delete gauss;
    }
    print::separator(10, ' ');
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
template <int D> void build_grid(FunctionTree<D> &out, FunctionTree<D> &inp, int maxIter) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, maxScale, nullptr);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
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
template <int D> void build_grid(FunctionTree<D> &out, FunctionTreeVector<D> &inp, int maxIter) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, maxScale, nullptr);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
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
template <int D> void copy_func(FunctionTree<D> &out, FunctionTree<D> &inp) {
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
template <int D> void copy_grid(FunctionTree<D> &out, FunctionTree<D> &inp) {
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
template <int D> void clear_grid(FunctionTree<D> &out) {
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
template <int D> int refine_grid(FunctionTree<D> &out, int scales) {
    auto nSplit = 0;
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    SplitAdaptor<D> adaptor(maxScale, true); // Splits all nodes
    for (auto n = 0; n < scales; n++) {
        nSplit += builder.split(out, adaptor, true); // Transfers coefs to children
    }
    return nSplit;
}

/** @brief Refine the grid of a MW function representation
 *
 * @param[in,out] out Output tree to be refined
 * @param[in] prec Precision for initial split check
 * @param[in] absPrec Build output tree based on absolute precision
 *
 * This will first perform a split check on the existing end nodes in the tree
 * based on the provided precision parameter, then it will compute scaling coefs
 * of the new nodes, thus leaving the function representation unchanged, but on a
 * larger grid.
 *
 */
template <int D> int refine_grid(FunctionTree<D> &out, double prec, bool absPrec) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale, absPrec);
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
template <int D> int refine_grid(FunctionTree<D> &out, FunctionTree<D> &inp) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, maxScale, nullptr);
    auto nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

template void build_grid(FunctionTree<1> &out, const Gaussian<1> &inp, int maxIter);
template void build_grid(FunctionTree<2> &out, const Gaussian<2> &inp, int maxIter);
template void build_grid(FunctionTree<3> &out, const Gaussian<3> &inp, int maxIter);

template void build_grid(FunctionTree<1> &out, const GaussExp<1> &inp, int maxIter);
template void build_grid(FunctionTree<2> &out, const GaussExp<2> &inp, int maxIter);
template void build_grid(FunctionTree<3> &out, const GaussExp<3> &inp, int maxIter);
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
template int refine_grid(FunctionTree<1> &out, double prec, bool absPrec);
template int refine_grid(FunctionTree<2> &out, double prec, bool absPrec);
template int refine_grid(FunctionTree<3> &out, double prec, bool absPrec);
template int refine_grid(FunctionTree<1> &out, FunctionTree<1> &inp);
template int refine_grid(FunctionTree<2> &out, FunctionTree<2> &inp);
template int refine_grid(FunctionTree<3> &out, FunctionTree<3> &inp);

} // namespace mrcpp
