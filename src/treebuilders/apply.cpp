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

template <int D, typename T> void apply_on_unit_cell(bool inside, double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, int maxIter, bool absPrec);

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
template <int D, typename T> void apply(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, int maxIter, bool absPrec) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");

    Timer pre_t;
    oper.calcBandWidths(prec);
    int maxScale = out.getMRA().getMaxScale();
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    ConvolutionCalculator<D, T> calculator(prec, oper, inp);
    pre_t.stop();
    TreeBuilder<D, T> builder;
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


/** @brief Application of MW integral convolution operator on Four component
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] oper: Convolution operator to apply
 * @param[in] inp: Input function
 * @param[in] metric: 4x4 array with coefficients that relates the in and out components
 * @param[in] maxIter: Maximum number of refinement iterations in output tree, default -1
 * @param[in] absPrec: Build output tree based on absolute precision, default false
 *
 * @details The output function will be computed using the general algorithm:
 * - For each input component apply the operator
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 * - After application multiply by metric coefficient, and put in relevant output component
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D, typename T> void apply(double prec, CompFunction<D, T> &out, ConvolutionOperator<D> &oper, CompFunction<D, T> &inp, T **metric, int maxIter, bool absPrec) {
    for (int icomp = 0; icomp < 4; icomp++){
        if (inp.Comp[icomp]!=nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++){
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    apply(prec, *out.Comp[ocomp], oper, *inp.Comp[icomp], maxIter, absPrec);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.Comp[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
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
template <int D, typename T> void apply_on_unit_cell(bool inside, double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, int maxIter, bool absPrec) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");

    Timer pre_t;
    oper.calcBandWidths(prec);
    int maxScale = out.getMRA().getMaxScale();
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    ConvolutionCalculator<D, T> calculator(prec, oper, inp);
    calculator.startManipulateOperator(inside);
    pre_t.stop();

    TreeBuilder<D, T> builder;
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
template <int D, typename T> void apply(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, FunctionTreeVector<D, T> &precTrees, int maxIter, bool absPrec) {
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

    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    adaptor.setPrecFunction(precFunc);
    ConvolutionCalculator<D, T> calculator(prec, oper, inp);
    calculator.setPrecFunction(precFunc);
    pre_t.stop();

    TreeBuilder<D, T> builder;
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

template <int D, typename T> void apply(double prec, CompFunction<D, T> &out, ConvolutionOperator<D> &oper, CompFunction<D, T> &inp, FunctionTreeVector<D, T> *precTrees, T **metric, int maxIter, bool absPrec) {
    for (int icomp = 0; icomp < 4; icomp++){
        if (inp.Comp[icomp]!=nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++){
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    apply(prec, *out.Comp[ocomp], oper, *inp.Comp[icomp], precTrees[icomp], maxIter, absPrec);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.Comp[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
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
template <int D, typename T> void apply_far_field(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, int maxIter, bool absPrec) {
    apply_on_unit_cell<D>(false, prec, out, oper, inp, maxIter, absPrec);
}

template <int D, typename T> void apply_far_field(double prec, CompFunction<D, T> &out, ConvolutionOperator<D> &oper, CompFunction<D, T> &inp, T **metric, int maxIter, bool absPrec) {
    for (int icomp = 0; icomp < 4; icomp++){
        if (inp.Comp[icomp]!=nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++){
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    apply_on_unit_cell<D>(false, prec, *out.Comp[ocomp], oper, *inp.Comp[icomp], maxIter, absPrec);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.Comp[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
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
template <int D, typename T> void apply_near_field(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper, FunctionTree<D, T> &inp, int maxIter, bool absPrec) {
    apply_on_unit_cell<D>(true, prec, out, oper, inp, maxIter, absPrec);
}


template <int D, typename T> void apply_near_field(double prec, CompFunction<D, T> &out, ConvolutionOperator<D> &oper, CompFunction<D, T> &inp, T **metric, int maxIter, bool absPrec) {
    for (int icomp = 0; icomp < 4; icomp++){
        if (inp.Comp[icomp]!=nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++){
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    apply_on_unit_cell<D>(true, prec, *out.Comp[ocomp], oper, *inp.Comp[icomp], maxIter, absPrec);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.Comp[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
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
template <int D, typename T> void apply(FunctionTree<D, T> &out, DerivativeOperator<D> &oper, FunctionTree<D, T> &inp, int dir) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    TreeBuilder<D, T> builder;
    int maxScale = out.getMRA().getMaxScale();

    int bw[D]; // Operator bandwidth in [x,y,z]
    for (int d = 0; d < D; d++) bw[d] = 0;

    // Copy input tree plus bandwidth in operator direction
    Timer pre_t;
    oper.calcBandWidths(1.0); // Fixed 0 or 1 for derivatives
    bw[dir] = oper.getMaxBandWidth();
    CopyAdaptor<D, T> pre_adaptor(inp, maxScale, bw);
    DefaultCalculator<D, T> pre_calculator;
    builder.build(out, pre_calculator, pre_adaptor, -1);
    pre_t.stop();

    // Apply operator on fixed expanded grid
    SplitAdaptor<D, T> apply_adaptor(maxScale, false); // Splits no nodes
    DerivativeCalculator<D, T> apply_calculator(dir, oper, inp);
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

template <int D, typename T> void apply(CompFunction<D, T> &out, DerivativeOperator<D> &oper, CompFunction<D, T> &inp, T **metric, int dir) {
    for (int icomp = 0; icomp < 4; icomp++){
        if (inp.Comp[icomp]!=nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++){
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    apply(*out.Comp[ocomp], oper, *inp.Comp[icomp], dir);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.Comp[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
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
template <int D, typename T> FunctionTreeVector<D, T> gradient(DerivativeOperator<D> &oper, FunctionTree<D, T> &inp) {
    FunctionTreeVector<D, T> out;
    for (int d = 0; d < D; d++) {
        auto *grad_d = new FunctionTree<D, T>(inp.getMRA());
        apply(*grad_d, oper, inp, d);
        out.push_back({1.0, grad_d});
    }
    return out;
}

template <int D, typename T> CompFunctionVector<D, T> gradient(DerivativeOperator<D> &oper, CompFunction<D, T> &inp, T **metric) {
    CompFunctionVector<D, T> out;
    for (int d = 0; d < D; d++) {
        CompFunction<D, T> *grad_d = new CompFunction<D, T>();
        for (int icomp = 0; icomp < 4; icomp++){
            if (inp.Comp[icomp]!=nullptr) {
                for (int ocomp = 0; ocomp < 4; ocomp++){
                    if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                        grad_d->Comp[ocomp] = new FunctionTree<D, T>(inp.getMRA());
                        apply(grad_d->Comp[ocomp], oper, *inp.Comp[icomp], d);
                        if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                            grad_d->Comp[ocomp]->rescale(metric[icomp][ocomp]);
                        }
                    }
                }
            }
        }
        out.oush_back(grad_d);
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
template <int D, typename T> void divergence(FunctionTree<D, T> &out, DerivativeOperator<D> &oper, FunctionTreeVector<D, T> &inp) {
    if (inp.size() != D) MSG_ABORT("Dimension mismatch");
    for (auto i = 0; i < inp.size(); i++)
        if (out.getMRA() != get_func(inp, i).getMRA()) MSG_ABORT("Incompatible MRA");

    FunctionTreeVector<D, T> tmp_vec;
    for (int d = 0; d < D; d++) {
        double coef_d = get_coef(inp, d);
        FunctionTree<D, T> &func_d = get_func(inp, d);
        auto *out_d = new FunctionTree<D, T>(func_d.getMRA());
        apply(*out_d, oper, func_d, d);
        tmp_vec.push_back(std::make_tuple(coef_d, out_d));
    }
    build_grid(out, tmp_vec);
    add(-1.0, out, tmp_vec, 0); // Addition on union grid
    clear(tmp_vec, true);
}

template <int D, typename T> void divergence(CompFunction<D, T> &out, DerivativeOperator<D> &oper, FunctionTreeVector<D, T> *inp, T **metric) {
    for (int icomp = 0; icomp < 4; icomp++){
        if (inp[icomp]!=nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++){
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    divergence(*out.Comp[ocomp], oper, inp[icomp]);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.Comp[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
}

template <int D, typename T> void divergence(FunctionTree<D, T> &out, DerivativeOperator<D> &oper, std::vector<FunctionTree<D, T> *> &inp) {
    FunctionTreeVector<D, T> inp_vec;
    for (auto &t : inp) inp_vec.push_back({1.0, t});
    divergence(out, oper, inp_vec);
}
template <int D, typename T> void divergence(CompFunction<D, T> &out, DerivativeOperator<D> &oper, std::vector<FunctionTree<D, T> *> *inp, T **metric) {
    for (int icomp = 0; icomp < 4; icomp++){
        if (inp[icomp]!=nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++){
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    apply(*out.Comp[ocomp], oper, inp[icomp]);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.Comp[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
}

template void apply<1, double>(double prec, FunctionTree<1, double> &out, ConvolutionOperator<1> &oper, FunctionTree<1, double> &inp, int maxIter, bool absPrec);
template void apply<2, double>(double prec, FunctionTree<2, double> &out, ConvolutionOperator<2> &oper, FunctionTree<2, double> &inp, int maxIter, bool absPrec);
template void apply<3, double>(double prec, FunctionTree<3, double> &out, ConvolutionOperator<3> &oper, FunctionTree<3, double> &inp, int maxIter, bool absPrec);
template void apply<1, double>(double prec, CompFunction<1, double> &out, ConvolutionOperator<1> &oper, CompFunction<1, double> &inp, double **metric, int maxIter = -1, bool absPrec = false);
template void apply<2, double>(double prec, CompFunction<2, double> &out, ConvolutionOperator<2> &oper, CompFunction<2, double> &inp, double **metric, int maxIter = -1, bool absPrec = false);
template void apply<3, double>(double prec, CompFunction<3, double> &out, ConvolutionOperator<3> &oper, CompFunction<3, double> &inp, double **metric, int maxIter = -1, bool absPrec = false);
template void apply<1, double>(double prec, FunctionTree<1, double> &out, ConvolutionOperator<1> &oper, FunctionTree<1, double> &inp, FunctionTreeVector<1, double> &precTrees, int maxIter, bool absPrec);
template void apply<2, double>(double prec, FunctionTree<2, double> &out, ConvolutionOperator<2> &oper, FunctionTree<2, double> &inp, FunctionTreeVector<2, double> &precTrees, int maxIter, bool absPrec);
template void apply<3, double>(double prec, FunctionTree<3, double> &out, ConvolutionOperator<3> &oper, FunctionTree<3, double> &inp, FunctionTreeVector<3, double> &precTrees, int maxIter, bool absPrec);
template void apply_far_field<1, double>(double prec, FunctionTree<1, double> &out, ConvolutionOperator<1> &oper, FunctionTree<1, double> &inp, int maxIter, bool absPrec);
template void apply_far_field<2, double>(double prec, FunctionTree<2, double> &out, ConvolutionOperator<2> &oper, FunctionTree<2, double> &inp, int maxIter, bool absPrec);
template void apply_far_field<3, double>(double prec, FunctionTree<3, double> &out, ConvolutionOperator<3> &oper, FunctionTree<3, double> &inp, int maxIter, bool absPrec);
template void apply_near_field<1, double>(double prec, FunctionTree<1, double> &out, ConvolutionOperator<1> &oper, FunctionTree<1, double> &inp, int maxIter, bool absPrec);
template void apply_near_field<2, double>(double prec, FunctionTree<2, double> &out, ConvolutionOperator<2> &oper, FunctionTree<2, double> &inp, int maxIter, bool absPrec);
template void apply_near_field<3, double>(double prec, FunctionTree<3, double> &out, ConvolutionOperator<3> &oper, FunctionTree<3, double> &inp, int maxIter, bool absPrec);
template void apply<1, double>(FunctionTree<1, double> &out, DerivativeOperator<1> &oper, FunctionTree<1, double> &inp, int dir);
template void apply<2, double>(FunctionTree<2, double> &out, DerivativeOperator<2> &oper, FunctionTree<2, double> &inp, int dir);
template void apply<3, double>(FunctionTree<3, double> &out, DerivativeOperator<3> &oper, FunctionTree<3, double> &inp, int dir);
template void divergence<1, double>(FunctionTree<1, double> &out, DerivativeOperator<1> &oper, FunctionTreeVector<1, double> &inp);
template void divergence<2, double>(FunctionTree<2, double> &out, DerivativeOperator<2> &oper, FunctionTreeVector<2, double> &inp);
template void divergence<3, double>(FunctionTree<3, double> &out, DerivativeOperator<3> &oper, FunctionTreeVector<3, double> &inp);
template void divergence<1, double>(FunctionTree<1, double> &out, DerivativeOperator<1> &oper, std::vector<FunctionTree<1, double> *> &inp);
template void divergence<2, double>(FunctionTree<2, double> &out, DerivativeOperator<2> &oper, std::vector<FunctionTree<2, double> *> &inp);
template void divergence<3, double>(FunctionTree<3, double> &out, DerivativeOperator<3> &oper, std::vector<FunctionTree<3, double> *> &inp);
template FunctionTreeVector<1, double> gradient<1>(DerivativeOperator<1> &oper, FunctionTree<1, double> &inp);
template FunctionTreeVector<2, double> gradient<2>(DerivativeOperator<2> &oper, FunctionTree<2, double> &inp);
template FunctionTreeVector<3, double> gradient<3>(DerivativeOperator<3> &oper, FunctionTree<3, double> &inp);



template void apply<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, ConvolutionOperator<1> &oper, FunctionTree<1, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, ConvolutionOperator<2> &oper, FunctionTree<2, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, ConvolutionOperator<3> &oper, FunctionTree<3, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply<1, ComplexDouble>(double prec, CompFunction<1, ComplexDouble> &out, ConvolutionOperator<1> &oper, CompFunction<1, ComplexDouble> &inp, ComplexDouble **metric, int maxIter = -1, bool absPrec = false);
template void apply<2, ComplexDouble>(double prec, CompFunction<2, ComplexDouble> &out, ConvolutionOperator<2> &oper, CompFunction<2, ComplexDouble> &inp, ComplexDouble **metric, int maxIter = -1, bool absPrec = false);
template void apply<3, ComplexDouble>(double prec, CompFunction<3, ComplexDouble> &out, ConvolutionOperator<3> &oper, CompFunction<3, ComplexDouble> &inp, ComplexDouble **metric, int maxIter = -1, bool absPrec = false);
template void apply<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, ConvolutionOperator<1> &oper, FunctionTree<1, ComplexDouble> &inp, FunctionTreeVector<1, ComplexDouble> &precTrees, int maxIter, bool absPrec);
template void apply<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, ConvolutionOperator<2> &oper, FunctionTree<2, ComplexDouble> &inp, FunctionTreeVector<2, ComplexDouble> &precTrees, int maxIter, bool absPrec);
template void apply<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, ConvolutionOperator<3> &oper, FunctionTree<3, ComplexDouble> &inp, FunctionTreeVector<3, ComplexDouble> &precTrees, int maxIter, bool absPrec);
template void apply_far_field<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, ConvolutionOperator<1> &oper, FunctionTree<1, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply_far_field<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, ConvolutionOperator<2> &oper, FunctionTree<2, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply_far_field<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, ConvolutionOperator<3> &oper, FunctionTree<3, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply_near_field<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, ConvolutionOperator<1> &oper, FunctionTree<1, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply_near_field<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, ConvolutionOperator<2> &oper, FunctionTree<2, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply_near_field<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, ConvolutionOperator<3> &oper, FunctionTree<3, ComplexDouble> &inp, int maxIter, bool absPrec);
template void apply<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, DerivativeOperator<1> &oper, FunctionTree<1, ComplexDouble> &inp, int dir);
template void apply<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, DerivativeOperator<2> &oper, FunctionTree<2, ComplexDouble> &inp, int dir);
template void apply<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, DerivativeOperator<3> &oper, FunctionTree<3, ComplexDouble> &inp, int dir);
template void divergence<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, DerivativeOperator<1> &oper, FunctionTreeVector<1, ComplexDouble> &inp);
template void divergence<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, DerivativeOperator<2> &oper, FunctionTreeVector<2, ComplexDouble> &inp);
template void divergence<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, DerivativeOperator<3> &oper, FunctionTreeVector<3, ComplexDouble> &inp);
template void divergence<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &out, DerivativeOperator<1> &oper, std::vector<FunctionTree<1, ComplexDouble> *> &inp);
template void divergence<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &out, DerivativeOperator<2> &oper, std::vector<FunctionTree<2, ComplexDouble> *> &inp);
template void divergence<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &out, DerivativeOperator<3> &oper, std::vector<FunctionTree<3, ComplexDouble> *> &inp);
template FunctionTreeVector<1, ComplexDouble> gradient<1>(DerivativeOperator<1> &oper, FunctionTree<1, ComplexDouble> &inp);
template FunctionTreeVector<2, ComplexDouble> gradient<2>(DerivativeOperator<2> &oper, FunctionTree<2, ComplexDouble> &inp);
template FunctionTreeVector<3, ComplexDouble> gradient<3>(DerivativeOperator<3> &oper, FunctionTree<3, ComplexDouble> &inp);

} // namespace mrcpp
