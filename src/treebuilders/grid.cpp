/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
#include "functions/Gaussian.h"
#include "functions/function_utils.h"
#include "utils/Printer.h"

namespace mrcpp {

/** @brief Build empty grid by uniform refinement
 *
 * @param[in,out] out: Output tree to be built
 * @param[in] scales: Number of refinement levels
 *
 * @details This will split ALL leaf nodes in the tree the given number of times.
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called.
 */
template <int D, typename T> void build_grid(FunctionTree<D, T> &out, int scales) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    DefaultCalculator<D, T> calculator;
    SplitAdaptor<D, T> adaptor(maxScale, true); // Splits all nodes
    for (auto n = 0; n < scales; n++) builder.build(out, calculator, adaptor, 1);
}

/** @brief Build empty grid based on info from analytic function
 *
 * @param[out] out: Output tree to be built
 * @param[in] inp: Input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 *
 * @details The grid of the output function will be EXTENDED using the general
 * algorithm:
 * - Loop through current leaf nodes of the output tree
 * - Refine node based on custom split check from the function
 * - Repeat until convergence or `maxIter` is reached
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called. It requires that the functions
 * `isVisibleAtScale()` and `isZeroOnInterval()` is implemented in the
 * particular `RepresentableFunction`.
 *
 */
template <int D, typename T> void build_grid(FunctionTree<D, T> &out, const RepresentableFunction<D, T> &inp, int maxIter) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    AnalyticAdaptor<D, T> adaptor(inp, maxScale);
    DefaultCalculator<D, T> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
}

/** @brief Build empty grid based on info from Gaussian expansion
 *
 * @param[out] out: Output tree to be built
 * @param[in] inp: Input Gaussian expansion
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 *
 * @details The grid of the output function will be EXTENDED using the general
 * algorithm:
 * - Loop through current leaf nodes of the output tree
 * - Refine node based on custom split check from the function
 * - Repeat until convergence or `maxIter` is reached
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called. It will loop through the Gaussians in the
 * expansion and extend the grid based on the position and exponent of each
 * term. Higher exponent means finer resolution.
 *
 */
template <int D> void build_grid(FunctionTree<D> &out, const GaussExp<D> &inp, int maxIter) {
    if (!out.getMRA().getWorldBox().isPeriodic()) {
        auto maxScale = out.getMRA().getMaxScale();
        TreeBuilder<D> builder;
        DefaultCalculator<D> calculator;
        for (auto i = 0; i < inp.size(); i++) {
            AnalyticAdaptor<D> adaptor(inp.getFunc(i), maxScale);
            builder.build(out, calculator, adaptor, maxIter);
        }
    } else {
        for (auto i = 0; i < inp.size(); i++) {
            auto *gauss = inp.getFunc(i).copy();
            build_grid(out, *gauss, maxIter);
            delete gauss;
        }
    }
    print::separator(10, ' ');
}

/** @brief Build empty grid based on another MW function representation
 *
 * @param[out] out: Output tree to be built
 * @param[in] inp: Input tree
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 *
 * @details The grid of the output function will be EXTENDED with all existing
 * nodes in corresponding input function, using the general algorithm:
 * - Loop through current leaf nodes of the output tree
 * - Refine node if the corresponding node in the input has children
 * - Repeat until all input nodes are covered or `maxIter` is reached
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called. This means that all nodes on the input
 * tree will also be in the final output tree (unless `maxIter` is reached,
 * but NOT vice versa.
 *
 */
template <int D, typename T> void build_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp, int maxIter) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    CopyAdaptor<D, T> adaptor(inp, maxScale, nullptr);
    DefaultCalculator<D, T> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
}

/** @brief Build empty grid based on several MW function representation
 *
 * @param[out] out: Output tree to be built
 * @param[in] inp: Input tree vector
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 *
 * @details The grid of the output function will be EXTENDED with all existing
 * nodes in all corresponding input functions, using the general algorithm:
 * - Loop through current leaf nodes of the output tree
 * - Refine node if the corresponding node in one of the inputs has children
 * - Repeat until all input nodes are covered or `maxIter` is reached
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called. This means that the final output grid
 * will contain (at least) the union of the nodes of all input trees (unless
 * `maxIter` is reached).
 *
 */
template <int D, typename T> void build_grid(FunctionTree<D, T> &out, FunctionTreeVector<D, T> &inp, int maxIter) {
    for (auto i = 0; i < inp.size(); i++)
        if (out.getMRA() != get_func(inp, i).getMRA()) MSG_ABORT("Incompatible MRA");

    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    CopyAdaptor<D, T> adaptor(inp, maxScale, nullptr);
    DefaultCalculator<D, T> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
}

template <int D, typename T> void build_grid(FunctionTree<D, T> &out, std::vector<FunctionTree<D, T> *> &inp, int maxIter) {
    FunctionTreeVector<D, T> inp_vec;
    for (auto *t : inp) inp_vec.push_back({1.0, t});
    build_grid(out, inp_vec, maxIter);
}

/** @brief Copy function from one tree onto the grid of another tree, fixed grid
 *
 * @param[out] out: Output function
 * @param[in] inp: Input function
 *
 * @details The output function will be computed using the general algorithm:
 * - Loop through current leaf nodes of the output tree
 * - Copy MW coefs from the corresponding input node
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called and will overwrite any existing coefs.
 *
 */
template <int D, typename T> void copy_func(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    FunctionTreeVector<D, T> tmp_vec;
    tmp_vec.push_back(std::make_tuple(1.0, &inp));
    add(-1.0, out, tmp_vec);
}

/** @brief Build empty grid that is identical to another MW grid
 *
 * @param[out] out: Output tree to be built
 * @param[in] inp: Input tree
 *
 * @note The difference from the corresponding `build_grid` function is that
 * this will first clear the grid of the `out` function, while `build_grid`
 * will _extend_ the existing grid.
 *
 */
template <int D, typename T> void copy_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA")
    out.clear();
    build_grid(out, inp);
}

/** @brief Build empty grid that is identical to another MW grid for every component
 *
 * @param[out] out: Output to be built
 * @param[in] inp: Input
 *
 * @note The difference from the corresponding `build_grid` function is that
 * this will first clear the grid of the `out` function, while `build_grid`
 * will _extend_ the existing grid.
 *
 */
template <int D> void copy_grid(CompFunction<D> &out, CompFunction<D> &inp) {
    out.free();
    out.func_ptr->data = inp.func_ptr->data;
    out.alloc(inp.Ncomp());
    for (int i = 0; i < inp.Ncomp(); i++) {
        if (inp.isreal()) build_grid(*out.CompD[i], *inp.CompD[i]);
        if (inp.iscomplex()) build_grid(*out.CompC[i], *inp.CompC[i]);
    }
}

/** @brief Clear the MW coefficients of a function representation
 *
 * @param[in,out] out: Output function to be cleared
 *
 * @note This will only clear the MW coefs in the existing nodes, it will not
 * change the grid of the function. Use `FunctionTree::clear()` to remove all
 * grid refinement as well.
 *
 */
template <int D, typename T> void clear_grid(FunctionTree<D, T> &out) {
    TreeBuilder<D, T> builder;
    DefaultCalculator<D, T> calculator;
    builder.clear(out, calculator);
}

/** @brief Refine the grid of a MW function representation
 *
 * @param[in,out] out: Output tree to be refined
 * @param[in] scales: Number of refinement levels
 * @returns The number of nodes that were split
 *
 * @details This will split ALL leaf nodes in the tree the given number of
 * times, then it will compute scaling coefs of the new nodes, thus leaving
 * the function representation unchanged, but on a larger grid.
 *
 */
template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, int scales) {
    auto nSplit = 0;
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    SplitAdaptor<D, T> adaptor(maxScale, true); // Splits all nodes
    for (auto n = 0; n < scales; n++) {
        nSplit += builder.split(out, adaptor, true); // Transfers coefs to children
    }
    return nSplit;
}

/** @brief Refine the grid of a MW function representation
 *
 * @param[in,out] out: Output tree to be refined
 * @param[in] prec: Precision for initial split check
 * @param[in] absPrec: Build output tree based on absolute precision
 * @returns The number of nodes that were split
 *
 * @details This will first perform a split check on the existing leaf nodes in
 * the tree based on the provided precision parameter, then it will compute
 * scaling coefs of the new nodes, thus leaving the function representation
 * unchanged, but (possibly) on a larger grid.
 *
 */
template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, double prec, bool absPrec) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    int nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

/** @brief Refine the grid of a MW function representation
 *
 * @param[in,out] out: Output tree to be refined
 * @param[in] inp: Input tree that defines the new grid
 * @returns The number of nodes that were split
 *
 * @details This will first perform a split check on the existing leaf nodes
 * in the output tree based on the structure of the input tree (same as
 * build_grid), then it will compute scaling coefs of the new nodes, thus
 * leaving the function representation unchanged, but on a larger grid.
 *
 */
template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA")
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    CopyAdaptor<D, T> adaptor(inp, maxScale, nullptr);
    auto nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

/** @brief Refine the grid of a MW function representation
 *
 * @param[in,out] out: Output tree to be refined
 * @param[in] inp: Input function
 *
 * @details This will first perform a split check on the existing leaf nodes
 * in the output tree based on the structure of the input function (same as
 * build_grid), then it will compute scaling coefs of the new nodes, thus
 * leaving the function representation unchanged, but on a larger grid.
 * It requires that the functions `isVisibleAtScale()` and `isZeroOnInterval()`
 * is implemented in the particular `RepresentableFunction`.
 *
 */
template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, const RepresentableFunction<D, T> &inp) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    AnalyticAdaptor<D, T> adaptor(inp, maxScale);
    int nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

template void copy_grid(CompFunction<1> &out, CompFunction<1> &inp);
template void copy_grid(CompFunction<2> &out, CompFunction<2> &inp);
template void copy_grid(CompFunction<3> &out, CompFunction<3> &inp);

template void build_grid<1, double>(FunctionTree<1, double> &out, int scales);
template void build_grid<2, double>(FunctionTree<2, double> &out, int scales);
template void build_grid<3, double>(FunctionTree<3, double> &out, int scales);
template void build_grid<1>(FunctionTree<1> &out, const GaussExp<1> &inp, int maxIter);
template void build_grid<2>(FunctionTree<2> &out, const GaussExp<2> &inp, int maxIter);
template void build_grid<3>(FunctionTree<3> &out, const GaussExp<3> &inp, int maxIter);
template void build_grid<1, double>(FunctionTree<1, double> &out, const RepresentableFunction<1, double> &inp, int maxIter);
template void build_grid<2, double>(FunctionTree<2, double> &out, const RepresentableFunction<2, double> &inp, int maxIter);
template void build_grid<3, double>(FunctionTree<3, double> &out, const RepresentableFunction<3, double> &inp, int maxIter);
template void build_grid<1, double>(FunctionTree<1, double> &out, FunctionTree<1, double> &inp, int maxIter);
template void build_grid<2, double>(FunctionTree<2, double> &out, FunctionTree<2, double> &inp, int maxIter);
template void build_grid<3, double>(FunctionTree<3, double> &out, FunctionTree<3, double> &inp, int maxIter);
template void build_grid<1, double>(FunctionTree<1, double> &out, FunctionTreeVector<1, double> &inp, int maxIter);
template void build_grid<2, double>(FunctionTree<2, double> &out, FunctionTreeVector<2, double> &inp, int maxIter);
template void build_grid<3, double>(FunctionTree<3, double> &out, FunctionTreeVector<3, double> &inp, int maxIter);
template void build_grid<1, double>(FunctionTree<1, double> &out, std::vector<FunctionTree<1, double> *> &inp, int maxIter);
template void build_grid<2, double>(FunctionTree<2, double> &out, std::vector<FunctionTree<2, double> *> &inp, int maxIter);
template void build_grid<3, double>(FunctionTree<3, double> &out, std::vector<FunctionTree<3, double> *> &inp, int maxIter);
template void copy_func<1, double>(FunctionTree<1, double> &out, FunctionTree<1, double> &inp);
template void copy_func<2, double>(FunctionTree<2, double> &out, FunctionTree<2, double> &inp);
template void copy_func<3, double>(FunctionTree<3, double> &out, FunctionTree<3, double> &inp);
template void copy_grid<1, double>(FunctionTree<1, double> &out, FunctionTree<1, double> &inp);
template void copy_grid<2, double>(FunctionTree<2, double> &out, FunctionTree<2, double> &inp);
template void copy_grid<3, double>(FunctionTree<3, double> &out, FunctionTree<3, double> &inp);
template void clear_grid<1, double>(FunctionTree<1, double> &out);
template void clear_grid<2, double>(FunctionTree<2, double> &out);
template void clear_grid<3, double>(FunctionTree<3, double> &out);
template int refine_grid<1, double>(FunctionTree<1, double> &out, int scales);
template int refine_grid<2, double>(FunctionTree<2, double> &out, int scales);
template int refine_grid<3, double>(FunctionTree<3, double> &out, int scales);
template int refine_grid<1, double>(FunctionTree<1, double> &out, double prec, bool absPrec);
template int refine_grid<2, double>(FunctionTree<2, double> &out, double prec, bool absPrec);
template int refine_grid<3, double>(FunctionTree<3, double> &out, double prec, bool absPrec);
template int refine_grid<1, double>(FunctionTree<1, double> &out, FunctionTree<1, double> &inp);
template int refine_grid<2, double>(FunctionTree<2, double> &out, FunctionTree<2, double> &inp);
template int refine_grid<3, double>(FunctionTree<3, double> &out, FunctionTree<3, double> &inp);
template int refine_grid<1, double>(FunctionTree<1, double> &out, const RepresentableFunction<1, double> &inp);
template int refine_grid<2, double>(FunctionTree<2, double> &out, const RepresentableFunction<2, double> &inp);
template int refine_grid<3, double>(FunctionTree<3, double> &out, const RepresentableFunction<3, double> &inp);

template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, int scales);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, int scales);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, int scales);
template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, const RepresentableFunction<1, ComplexDouble> &inp, int maxIter);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, const RepresentableFunction<2, ComplexDouble> &inp, int maxIter);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, const RepresentableFunction<3, ComplexDouble> &inp, int maxIter);
template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &inp, int maxIter);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &inp, int maxIter);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &inp, int maxIter);
template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTreeVector<1, ComplexDouble> &inp, int maxIter);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTreeVector<2, ComplexDouble> &inp, int maxIter);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTreeVector<3, ComplexDouble> &inp, int maxIter);
template void build_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, std::vector<FunctionTree<1, ComplexDouble> *> &inp, int maxIter);
template void build_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, std::vector<FunctionTree<2, ComplexDouble> *> &inp, int maxIter);
template void build_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, std::vector<FunctionTree<3, ComplexDouble> *> &inp, int maxIter);
template void copy_func<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &inp);
template void copy_func<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &inp);
template void copy_func<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &inp);
template void copy_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &inp);
template void copy_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &inp);
template void copy_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &inp);
template void clear_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out);
template void clear_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out);
template void clear_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out);
template int refine_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, int scales);
template int refine_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, int scales);
template int refine_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, int scales);
template int refine_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, double prec, bool absPrec);
template int refine_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, double prec, bool absPrec);
template int refine_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, double prec, bool absPrec);
template int refine_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &inp);
template int refine_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &inp);
template int refine_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &inp);
template int refine_grid<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, const RepresentableFunction<1, ComplexDouble> &inp);
template int refine_grid<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, const RepresentableFunction<2, ComplexDouble> &inp);
template int refine_grid<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, const RepresentableFunction<3, ComplexDouble> &inp);

} // namespace mrcpp
