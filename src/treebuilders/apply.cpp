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

#include "apply.h"
#include "ConvolutionCalculator.h"
#include "CopyAdaptor.h"
#include "DefaultCalculator.h"
#include "DerivativeCalculator.h"
#include "SplitAdaptor.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "grid.h"
#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

template <int D> void apply_on_unit_cell(bool inside, double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter, bool absPrec);

/** @brief Application of MW integral convolution operator
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] oper: Convolution operator to apply
 * @param[in] inp: Input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree, default -1
 * @param[in] absPrec: Build output tree based on absolute precision, default false
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D> void apply(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter, bool absPrec) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");

    Timer pre_t;
    oper.calcBandWidths(prec);
    int maxScale = out.getMRA().getMaxScale();
    WaveletAdaptor<D> adaptor(prec, maxScale, absPrec);
    ConvolutionCalculator<D> calculator(prec, oper, inp);
    pre_t.stop();
    TreeBuilder<D> builder;
    builder.build(out, calculator, adaptor, maxIter);

    Timer post_t;
    oper.clearBandWidths();
    out.mwTransform(TopDown, false); // add coarse scale contributions
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    out.deleteGeneratedParents();
    inp.deleteGenerated();
    inp.deleteGeneratedParents();
    post_t.stop();

    print::time(10, "Time pre operator", pre_t);
    print::time(10, "Time post operator", post_t);
    print::separator(10, ' ');
}

/** @brief Application of MW integral convolution operator
 *
 * @param[in] inside: Use points inside (true) or outside (false) the unitcell
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] oper: Convolution operator to apply
 * @param[in] inp: Input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree, default -1
 * @param[in] absPrec: Build output tree based on absolute precision, default false
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D> void apply_on_unit_cell(bool inside, double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter, bool absPrec) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");

    Timer pre_t;
    oper.calcBandWidths(prec);
    int maxScale = out.getMRA().getMaxScale();
    WaveletAdaptor<D> adaptor(prec, maxScale, absPrec);
    ConvolutionCalculator<D> calculator(prec, oper, inp);
    calculator.startManipulateOperator(inside);
    pre_t.stop();

    TreeBuilder<D> builder;
    builder.build(out, calculator, adaptor, maxIter);

    Timer post_t;
    oper.clearBandWidths();
    out.mwTransform(TopDown, false); // add coarse scale contributions
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    out.deleteGeneratedParents();
    inp.deleteGenerated();
    inp.deleteGeneratedParents();
    post_t.stop();

    print::time(10, "Time pre operator", pre_t);
    print::time(10, "Time post operator", post_t);
    print::separator(10, ' ');
}

/** @brief Application of MW integral convolution operator
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] oper: Convolution operator to apply
 * @param[in] inp: Input function
 * @param[in] precTrees: Precision trees
 * @param[in] maxIter: Maximum number of refinement iterations in output tree, default -1
 * @param[in] absPrec: Build output tree based on absolute precision, default false
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on _scaled_ `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 *
 * The precision will be scaled locally by the maxNorms of the precTrees input vector.
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D> void apply(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, FunctionTreeVector<D> &precTrees, int maxIter, bool absPrec) {
    Timer pre_t;
    oper.calcBandWidths(prec);
    int maxScale = out.getMRA().getMaxScale();

    // The local precision will be scaled by the maxNorm of the
    // corresponding node(s) in the precTrees vector.
    for (int i = 0; i < precTrees.size(); i++) get_func(precTrees, i).makeMaxSquareNorms();
    auto precFunc = [&precTrees](const NodeIndex<D> &idx) -> double {
        auto maxNorm = (precTrees.size()) ? 0.0 : 1.0;
        for (int i = 0; i < precTrees.size(); i++) {
            auto &pNode = get_func(precTrees, i).getNode(idx);
            maxNorm = std::max(maxNorm, std::sqrt(pNode.getMaxSquareNorm()));
        }
        return 1.0 / maxNorm;
    };

    WaveletAdaptor<D> adaptor(prec, maxScale, absPrec);
    adaptor.setPrecFunction(precFunc);
    ConvolutionCalculator<D> calculator(prec, oper, inp);
    calculator.setPrecFunction(precFunc);
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

    print::time(10, "Time pre operator", pre_t);
    print::time(10, "Time post operator", post_t);
    print::separator(10, ' ');
}

/** @brief Application of MW integral convolution operator on a periodic cell,
           excluding contributions inside the unit cell.
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] oper: Convolution operator to apply
 * @param[in] inp: Input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree, default -1
 * @param[in] absPrec: Build output tree based on absolute precision, default false
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D> void apply_far_field(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter, bool absPrec) {
    apply_on_unit_cell<D>(false, prec, out, oper, inp, maxIter, absPrec);
}

/** @brief Application of MW integral convolution operator on a periodic cell,
           excluding contributions outside the unit cell.
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] oper: Convolution operator to apply
 * @param[in] inp: Input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree, default -1
 * @param[in] absPrec: Build output tree based on absolute precision, default false
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D> void apply_near_field(double prec, FunctionTree<D> &out, ConvolutionOperator<D> &oper, FunctionTree<D> &inp, int maxIter, bool absPrec) {
    apply_on_unit_cell<D>(true, prec, out, oper, inp, maxIter, absPrec);
}

/** @brief Application of MW derivative operator
 *
 * @param[out] out: Output function to be built
 * @param[in] oper: Derivative operator to apply
 * @param[in] inp: Input function
 * @param[in] dir: Direction of derivative
 *
 * @details The output function will be computed on a FIXED grid that is
 * predetermined by the type of derivative operator. For a strictly local
 * operator (ABGV_00), the grid is an exact copy of the input function. For
 * operators that involve also neighboring nodes (ABGV_55, PH, BS) the base grid
 * will be WIDENED by one node in the direction of application (on each side).
 *
 * @note The output function should contain only empty root nodes at entry.
 *
 */
template <int D> void apply(FunctionTree<D> &out, DerivativeOperator<D> &oper, FunctionTree<D> &inp, int dir) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
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
    if (out.isPeriodic()) out.rescale(std::pow(2.0, -oper.getOperatorRoot()));

    Timer post_t;
    oper.clearBandWidths();
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    inp.deleteGenerated();
    post_t.stop();

    print::time(10, "Time pre operator", pre_t);
    print::time(10, "Time post operator", post_t);
    print::separator(10, ' ');
}

/** @brief Calculation of gradient vector of a function
 *
 * @param[in] oper: Derivative operator to apply
 * @param[in] inp: Input function
 * @returns FunctionTreeVector containing the gradient
 *
 * @details The derivative operator is applied in each Cartesian direction to
 * the input function and appended to the output vector.
 *
 * @note The length of the output vector will be the template dimension D.
 *
 */
template <int D> FunctionTreeVector<D> gradient(DerivativeOperator<D> &oper, FunctionTree<D> &inp) {
    FunctionTreeVector<D> out;
    for (int d = 0; d < D; d++) {
        auto *grad_d = new FunctionTree<D>(inp.getMRA());
        apply(*grad_d, oper, inp, d);
        out.push_back({1.0, grad_d});
    }
    return out;
}

/** @brief Calculation of divergence of a function vector
 *
 * @param[out] out: Output function
 * @param[in] oper: Derivative operator to apply
 * @param[in] inp: Input function vector
 *
 * @details The derivative operator is applied in each Cartesian direction to
 * the corresponding components of the input vector and added up to the final
 * output. The grid of the output is fixed as the union of the component
 * grids (including any derivative widening, see derivative apply).
 *
 * @note
 * - The length of the input vector must be the same as the template dimension D.
 * - The output function should contain only empty root nodes at entry.
 *
 */
template <int D> void divergence(FunctionTree<D> &out, DerivativeOperator<D> &oper, FunctionTreeVector<D> &inp) {
    if (inp.size() != D) MSG_ABORT("Dimension mismatch");
    for (auto i = 0; i < inp.size(); i++)
        if (out.getMRA() != get_func(inp, i).getMRA()) MSG_ABORT("Incompatible MRA");

    FunctionTreeVector<D> tmp_vec;
    for (int d = 0; d < D; d++) {
        double coef_d = get_coef(inp, d);
        FunctionTree<D> &func_d = get_func(inp, d);
        auto *out_d = new FunctionTree<D>(func_d.getMRA());
        apply(*out_d, oper, func_d, d);
        tmp_vec.push_back(std::make_tuple(coef_d, out_d));
    }
    build_grid(out, tmp_vec);
    add(-1.0, out, tmp_vec, 0); // Addition on union grid
    clear(tmp_vec, true);
}

template <int D> void divergence(FunctionTree<D> &out, DerivativeOperator<D> &oper, std::vector<FunctionTree<D> *> &inp) {
    FunctionTreeVector<D> inp_vec;
    for (auto &t : inp) inp_vec.push_back({1.0, t});
    divergence(out, oper, inp_vec);
}

template void apply<1>(double prec, FunctionTree<1> &out, ConvolutionOperator<1> &oper, FunctionTree<1> &inp, int maxIter, bool absPrec);
template void apply<2>(double prec, FunctionTree<2> &out, ConvolutionOperator<2> &oper, FunctionTree<2> &inp, int maxIter, bool absPrec);
template void apply<3>(double prec, FunctionTree<3> &out, ConvolutionOperator<3> &oper, FunctionTree<3> &inp, int maxIter, bool absPrec);
template void apply<1>(double prec, FunctionTree<1> &out, ConvolutionOperator<1> &oper, FunctionTree<1> &inp, FunctionTreeVector<1> &precTrees, int maxIter, bool absPrec);
template void apply<2>(double prec, FunctionTree<2> &out, ConvolutionOperator<2> &oper, FunctionTree<2> &inp, FunctionTreeVector<2> &precTrees, int maxIter, bool absPrec);
template void apply<3>(double prec, FunctionTree<3> &out, ConvolutionOperator<3> &oper, FunctionTree<3> &inp, FunctionTreeVector<3> &precTrees, int maxIter, bool absPrec);
template void apply_far_field<1>(double prec, FunctionTree<1> &out, ConvolutionOperator<1> &oper, FunctionTree<1> &inp, int maxIter, bool absPrec);
template void apply_far_field<2>(double prec, FunctionTree<2> &out, ConvolutionOperator<2> &oper, FunctionTree<2> &inp, int maxIter, bool absPrec);
template void apply_far_field<3>(double prec, FunctionTree<3> &out, ConvolutionOperator<3> &oper, FunctionTree<3> &inp, int maxIter, bool absPrec);
template void apply_near_field<1>(double prec, FunctionTree<1> &out, ConvolutionOperator<1> &oper, FunctionTree<1> &inp, int maxIter, bool absPrec);
template void apply_near_field<2>(double prec, FunctionTree<2> &out, ConvolutionOperator<2> &oper, FunctionTree<2> &inp, int maxIter, bool absPrec);
template void apply_near_field<3>(double prec, FunctionTree<3> &out, ConvolutionOperator<3> &oper, FunctionTree<3> &inp, int maxIter, bool absPrec);
template void apply<1>(FunctionTree<1> &out, DerivativeOperator<1> &oper, FunctionTree<1> &inp, int dir);
template void apply<2>(FunctionTree<2> &out, DerivativeOperator<2> &oper, FunctionTree<2> &inp, int dir);
template void apply<3>(FunctionTree<3> &out, DerivativeOperator<3> &oper, FunctionTree<3> &inp, int dir);
template void divergence<1>(FunctionTree<1> &out, DerivativeOperator<1> &oper, FunctionTreeVector<1> &inp);
template void divergence<2>(FunctionTree<2> &out, DerivativeOperator<2> &oper, FunctionTreeVector<2> &inp);
template void divergence<3>(FunctionTree<3> &out, DerivativeOperator<3> &oper, FunctionTreeVector<3> &inp);
template void divergence<1>(FunctionTree<1> &out, DerivativeOperator<1> &oper, std::vector<FunctionTree<1> *> &inp);
template void divergence<2>(FunctionTree<2> &out, DerivativeOperator<2> &oper, std::vector<FunctionTree<2> *> &inp);
template void divergence<3>(FunctionTree<3> &out, DerivativeOperator<3> &oper, std::vector<FunctionTree<3> *> &inp);
template FunctionTreeVector<1> gradient<1>(DerivativeOperator<1> &oper, FunctionTree<1> &inp);
template FunctionTreeVector<2> gradient<2>(DerivativeOperator<2> &oper, FunctionTree<2> &inp);
template FunctionTreeVector<3> gradient<3>(DerivativeOperator<3> &oper, FunctionTree<3> &inp);

} // namespace mrcpp
