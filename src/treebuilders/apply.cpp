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

/**
 * @file apply.cpp
 * @brief Application pipelines for MW operators (convolution and derivative) to MW function trees.
 *
 * @details
 * This module provides high-level procedures to **apply multiresolution operators**
 * to MW representations of functions. Two broad operator families are supported:
 *
 * - **Convolution-like integral operators** (e.g., Poisson, Helmholtz, Heat, identity),
 *   implemented as separable kernels in the scaling basis via @ref mrcpp::ConvolutionOperator.
 *   Application is performed on an **adaptively refined** output grid to meet a target precision.
 *
 * - **Local or band-limited derivative operators** (e.g., ABGV, PH, BS) implemented via
 *   @ref mrcpp::DerivativeOperator. Application occurs on a **fixed grid** derived from the
 *   input and widened according to the operator bandwidth in the selected direction.
 *
 * The typical adaptive application pipeline for convolution operators is:
 * - Pre-step: estimate operator bandwidths at each scale and set up an adaptive refinement policy.
 * - Build-step: evaluate local operator actions on the current grid, refine where needed until
 *   the precision target is reached (or a maximum number of iterations is met).
 * - Post-step: assemble and transform the output to the MW domain, compute norms, and clean any
 *   transient data generated on the inputs.
 *
 * Additional features:
 * - **Near-/Far-field splits** on periodic domains by including/excluding contributions
 *   from the unit cell.
 * - **Precision scaling** using auxiliary trees that modulate local tolerances based on
 *   maximum norms.
 * - **Multi-component support** through a 4×4 metric that mixes input/output components
 *   for relativistic-like workflows.
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
#include <optional>

namespace mrcpp {

/**
 * @brief Internal helper to apply a convolution operator while restricting contributions
 *        to inside or outside of the unit cell on periodic domains.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type (e.g., double or ComplexDouble).
 *
 * @param[in]  inside   If true, include only contributions from inside the unit cell;
 *                      if false, include only contributions from outside the unit cell.
 * @param[in]  prec     Target precision that drives adaptive refinement.
 * @param[out] out      Output function to be built. Should contain empty root nodes on entry.
 * @param[in]  oper     Convolution operator to apply.
 * @param[in]  inp      Input function.
 * @param[in]  maxIter  Maximum number of refinement iterations. Negative means unbounded.
 * @param[in]  absPrec  If true, treat `prec` as an absolute tolerance; otherwise relative.
 *
 * @details
 * Follows the standard adaptive pipeline for convolution operators, with the difference that
 * the calculator is instructed to selectively include unit-cell contributions according to
 * the `inside` flag.
 */
template <int D, typename T>
void apply_on_unit_cell(bool inside,
                        double prec,
                        FunctionTree<D, T> &out,
                        ConvolutionOperator<D> &oper,
                        FunctionTree<D, T> &inp,
                        int maxIter,
                        bool absPrec);

/**
 * @brief Apply a convolution-like integral operator on a single-component function (adaptive).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  prec     Target precision driving adaptive refinement.
 * @param[out] out      Output function tree. Must belong to the same MRA as `inp`.
 * @param[in]  oper     Convolution operator to apply.
 * @param[in]  inp      Input function tree.
 * @param[in]  maxIter  Maximum refinement iterations (negative for unbounded, zero disables refinement).
 * @param[in]  absPrec  If true, treat `prec` as absolute; otherwise relative.
 *
 * @details
 * Pipeline:
 * - Pre: compute operator bandwidths and create a @ref WaveletAdaptor with the given precision policy.
 * - Build: @ref TreeBuilder iteratively refines and evaluates the operator action.
 * - Post: transform to MW domain, compute squared norms, and clean generated structures.
 *
 * @note The output tree should initially contain only empty root nodes.
 * @throws Aborts if `out` and `inp` belong to different MRAs.
 */
template <int D, typename T>
void apply(double prec,
           FunctionTree<D, T> &out,
           ConvolutionOperator<D> &oper,
           FunctionTree<D, T> &inp,
           int maxIter,
           bool absPrec) {
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
    out.mwTransform(TopDown, false);
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

/**
 * @brief Apply a convolution operator to a 4-component function using a mixing metric.
 *
 * @tparam D Spatial dimension.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output multi-component function (structure copied from `inp`).
 * @param[in]  oper     Convolution operator to apply.
 * @param[in]  inp      Input multi-component function.
 * @param[in]  metric   4×4 coefficient array mapping input to output components.
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute-vs-relative precision flag.
 *
 * @details
 * For each input component `icomp`, the operator is applied and accumulated into each output
 * component `ocomp` with weight `metric[icomp][ocomp]`. Real and complex specializations are
 * handled, including rescaling of the result by the metric entries.
 */
template <int D>
void apply(double prec,
           CompFunction<D> &out,
           ConvolutionOperator<D> &oper,
           const CompFunction<D> &inp,
           const ComplexDouble (*metric)[4],
           int maxIter,
           bool absPrec) {

    out = inp.paramCopy(true);
    for (int icomp = 0; icomp < inp.Ncomp(); icomp++) {
        for (int ocomp = 0; ocomp < 4; ocomp++) {
            if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                if (inp.isreal()) {
                    if (out.CompD[ocomp] == nullptr) out.alloc_comp(ocomp);
                    apply(prec, *out.CompD[ocomp], oper, *inp.CompD[icomp], maxIter, absPrec);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.CompD[ocomp]->rescale(metric[icomp][ocomp].real());
                    }
                } else {
                    if (out.CompC[ocomp] == nullptr) out.alloc_comp(ocomp);
                    apply(prec, *out.CompC[ocomp], oper, *inp.CompC[icomp], maxIter, absPrec);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.CompC[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
}

/**
 * @brief Apply a convolution operator while selectively including or excluding unit-cell contributions.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  inside   Select inside (true) or outside (false) contributions in the unit cell.
 * @param[in]  prec     Target precision.
 * @param[out] out      Output tree.
 * @param[in]  oper     Convolution operator.
 * @param[in]  inp      Input tree.
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute-vs-relative precision flag.
 */
template <int D, typename T>
void apply_on_unit_cell(bool inside,
                        double prec,
                        FunctionTree<D, T> &out,
                        ConvolutionOperator<D> &oper,
                        FunctionTree<D, T> &inp,
                        int maxIter,
                        bool absPrec) {
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
    out.mwTransform(TopDown, false);
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

/**
 * @brief Apply a convolution operator with **locally scaled precision** from auxiliary trees.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  prec       Base precision target.
 * @param[out] out        Output function tree.
 * @param[in]  oper       Convolution operator.
 * @param[in]  inp        Input function tree.
 * @param[in]  precTrees  Vector of trees whose max norms modulate the local precision.
 * @param[in]  maxIter    Maximum refinement iterations.
 * @param[in]  absPrec    Absolute-vs-relative precision flag.
 *
 * @details
 * The local precision at node index `idx` is scaled by `1 / max_norm(idx)`, where `max_norm`
 * is taken across the supplied `precTrees`. This provides an error budget that adapts to
 * local magnitudes of reference fields.
 */
template <int D, typename T>
void apply(double prec,
           FunctionTree<D, T> &out,
           ConvolutionOperator<D> &oper,
           FunctionTree<D, T> &inp,
           FunctionTreeVector<D, T> &precTrees,
           int maxIter,
           bool absPrec) {
    Timer pre_t;
    oper.calcBandWidths(prec);
    int maxScale = out.getMRA().getMaxScale();

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
    out.mwTransform(TopDown, false);
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    inp.deleteGenerated();
    post_t.stop();

    print::time(10, "Time pre operator", pre_t);
    print::time(10, "Time post operator", post_t);
    print::separator(10, ' ');
}

/**
 * @brief Multi-component variant of the precision-scaled convolution application.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  prec       Base precision target.
 * @param[out] out        Output multi-component function (structure copied from `inp`).
 * @param[in]  oper       Convolution operator.
 * @param[in]  inp        Input multi-component function.
 * @param[in]  precTrees  Array (per input component) of precision trees used for scaling.
 * @param[in]  metric     4×4 mixing matrix.
 * @param[in]  maxIter    Maximum refinement iterations.
 * @param[in]  absPrec    Absolute-vs-relative precision flag.
 */
template <int D, typename T>
void apply(double prec,
           CompFunction<D> &out,
           ConvolutionOperator<D> &oper,
           CompFunction<D> &inp,
           FunctionTreeVector<D, T> *precTrees,
           const ComplexDouble (*metric)[4],
           int maxIter,
           bool absPrec) {

    out = inp.paramCopy(true);
    for (int icomp = 0; icomp < inp.Ncomp(); icomp++) {
        for (int ocomp = 0; ocomp < 4; ocomp++) {
            if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                if (inp.isreal()) {
                    apply(prec, *out.CompD[ocomp], oper, *inp.CompD[icomp], precTrees[icomp], maxIter, absPrec);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.CompD[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                } else {
                    apply(prec, *out.CompC[ocomp], oper, *inp.CompC[icomp], precTrees[icomp], maxIter, absPrec);
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.CompC[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
}

/**
 * @brief Apply a convolution operator while excluding inside-cell contributions (far-field).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output function.
 * @param[in]  oper     Convolution operator.
 * @param[in]  inp      Input function.
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute-vs-relative precision flag.
 */
template <int D, typename T>
void apply_far_field(double prec,
                     FunctionTree<D, T> &out,
                     ConvolutionOperator<D> &oper,
                     FunctionTree<D, T> &inp,
                     int maxIter,
                     bool absPrec) {
    apply_on_unit_cell<D>(false, prec, out, oper, inp, maxIter, absPrec);
}

/**
 * @brief Multi-component far-field application with mixing metric.
 *
 * @tparam D Spatial dimension.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output multi-component function.
 * @param[in]  oper     Convolution operator.
 * @param[in]  inp      Input multi-component function.
 * @param[in]  metric   4×4 mixing matrix.
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute-vs-relative precision flag.
 */
template <int D>
void apply_far_field(double prec,
                     CompFunction<D> &out,
                     ConvolutionOperator<D> &oper,
                     CompFunction<D> &inp,
                     const ComplexDouble (*metric)[4],
                     int maxIter,
                     bool absPrec) {

    out = inp.paramCopy(true);
    for (int icomp = 0; icomp < 4; icomp++) {
        if (inp.Comp[icomp] != nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++) {
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    if (inp.isreal()) {
                        apply_on_unit_cell<D>(false, prec, *out.CompD[ocomp], oper, *inp.CompD[icomp], maxIter, absPrec);
                        if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                            out.CompD[ocomp]->rescale(metric[icomp][ocomp]);
                        }
                    } else {
                        apply_on_unit_cell<D>(false, prec, *out.CompC[ocomp], oper, *inp.CompC[icomp], maxIter, absPrec);
                        if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                            out.CompC[ocomp]->rescale(metric[icomp][ocomp]);
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Apply a convolution operator while excluding outside-cell contributions (near-field).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output function.
 * @param[in]  oper     Convolution operator.
 * @param[in]  inp      Input function.
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute-vs-relative precision flag.
 */
template <int D, typename T>
void apply_near_field(double prec,
                      FunctionTree<D, T> &out,
                      ConvolutionOperator<D> &oper,
                      FunctionTree<D, T> &inp,
                      int maxIter,
                      bool absPrec) {
    apply_on_unit_cell<D>(true, prec, out, oper, inp, maxIter, absPrec);
}

/**
 * @brief Multi-component near-field application with mixing metric.
 *
 * @tparam D Spatial dimension.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output multi-component function.
 * @param[in]  oper     Convolution operator.
 * @param[in]  inp      Input multi-component function.
 * @param[in]  metric   4×4 mixing matrix.
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute-vs-relative precision flag.
 */
template <int D>
void apply_near_field(double prec,
                      CompFunction<D> &out,
                      ConvolutionOperator<D> &oper,
                      CompFunction<D> &inp,
                      const ComplexDouble (*metric)[4],
                      int maxIter,
                      bool absPrec) {

    for (int icomp = 0; icomp < 4; icomp++) {
        if (inp.Comp[icomp] != nullptr) {
            for (int ocomp = 0; ocomp < 4; ocomp++) {
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    if (inp.isreal()) {
                        apply_on_unit_cell<D>(true, prec, *out.CompD[ocomp], oper, *inp.CompD[icomp], maxIter, absPrec);
                        if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                            out.CompD[ocomp]->rescale(metric[icomp][ocomp]);
                        }
                    } else {
                        apply_on_unit_cell<D>(true, prec, *out.CompC[ocomp], oper, *inp.CompC[icomp], maxIter, absPrec);
                        if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                            out.CompC[ocomp]->rescale(metric[icomp][ocomp]);
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Apply a **derivative operator** on a fixed grid in the given direction.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[out] out  Output function. Should contain only empty root nodes on entry.
 * @param[in]  oper Derivative operator (defines bandwidth and assembly policy).
 * @param[in]  inp  Input function.
 * @param[in]  dir  Direction of application (0 for x, 1 for y, 2 for z).
 *
 * @details
 * The output grid is constructed by copying the input grid and **widening** it by the
 * operator bandwidth along the selected direction, if needed. Application then proceeds
 * on this fixed grid without additional refinement.
 */
template <int D, typename T>
void apply(FunctionTree<D, T> &out,
           DerivativeOperator<D> &oper,
           FunctionTree<D, T> &inp,
           int dir) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    TreeBuilder<D, T> builder;
    int maxScale = out.getMRA().getMaxScale();

    int bw[D];
    for (int d = 0; d < D; d++) bw[d] = 0;

    Timer pre_t;
    oper.calcBandWidths(1.0);
    bw[dir] = oper.getMaxBandWidth();
    CopyAdaptor<D, T> pre_adaptor(inp, maxScale, bw);
    DefaultCalculator<D, T> pre_calculator;
    builder.build(out, pre_calculator, pre_adaptor, -1);
    pre_t.stop();

    SplitAdaptor<D, T> apply_adaptor(maxScale, false);
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

/**
 * @brief Multi-component derivative application with mixing metric.
 *
 * @tparam D Spatial dimension.
 *
 * @param[out] out    Output multi-component function.
 * @param[in]  oper   Derivative operator.
 * @param[in]  inp    Input multi-component function.
 * @param[in]  dir    Direction of derivative.
 * @param[in]  metric 4×4 mixing matrix.
 *
 * @details
 * Applies the derivative in `dir` to each input component and accumulates the result into
 * output components according to `metric`. Handles real-to-complex promotion if necessary.
 */
template <int D>
void apply(CompFunction<D> &out,
           DerivativeOperator<D> &oper,
           CompFunction<D> &inp,
           int dir,
           const ComplexDouble (*metric)[4]) {
    out = inp.paramCopy(true);
    for (int icomp = 0; icomp < inp.Ncomp(); icomp++) {
        for (int ocomp = 0; ocomp < 4; ocomp++) {
            if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                if (inp.isreal() && (std::imag(metric[icomp][ocomp]) < MachinePrec || inp.Ncomp() == 1)) {
                    apply(*out.CompD[ocomp], oper, *inp.CompD[icomp], dir);
                    if (std::norm(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        if (std::imag(metric[icomp][ocomp]) < MachinePrec)
                            out.CompD[ocomp]->rescale(std::real(metric[icomp][ocomp]));
                        else
                            out.func_ptr->data.c1[ocomp] *= metric[icomp][ocomp];
                    }
                    out.func_ptr->isreal = 1;
                } else {
                    if (inp.isreal()) {
                        apply(*out.CompD[ocomp], oper, *inp.CompD[icomp], dir);
                        out.CompD[icomp]->CopyTreeToComplex(out.CompC[ocomp]);
                        out.func_ptr->isreal = 0;
                        out.func_ptr->iscomplex = 1;
                    } else {
                        apply(*out.CompC[ocomp], oper, *inp.CompC[icomp], dir);
                    }
                    if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                        out.CompC[ocomp]->rescale(metric[icomp][ocomp]);
                    }
                }
            }
        }
    }
}

/**
 * @brief Compute the gradient vector of a scalar function using a derivative operator.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  oper Derivative operator to apply in each Cartesian direction.
 * @param[in]  inp  Input scalar function.
 * @return     FunctionTreeVector containing D components of the gradient.
 *
 * @details
 * Applies the operator in each direction `d = 0..D-1` and returns the resulting
 * component trees with unit weights.
 */
template <int D, typename T>
FunctionTreeVector<D, T> gradient(DerivativeOperator<D> &oper,
                                  FunctionTree<D, T> &inp) {
    FunctionTreeVector<D, T> out;
    for (int d = 0; d < D; d++) {
        auto *grad_d = new FunctionTree<D, T>(inp.getMRA());
        apply(*grad_d, oper, inp, d);
        out.push_back({1.0, grad_d});
    }
    return out;
}

/**
 * @brief Compute the gradient for 3D multi-component inputs with mixing metric.
 *
 * @param[in]  oper   Derivative operator.
 * @param[in]  inp    Input multi-component function.
 * @param[in]  metric 4×4 mixing matrix.
 * @return     Vector of component functions for each spatial direction.
 *
 * @details
 * For each spatial direction, applies the derivative operator to each component and
 * mixes according to `metric`. Handles both real and complex cases.
 */
std::vector<CompFunction<3> *> gradient(DerivativeOperator<3> &oper,
                                        CompFunction<3> &inp,
                                        const ComplexDouble (*metric)[4]) {
    std::vector<CompFunction<3> *> out;

    for (int d = 0; d < 3; d++) {
        CompFunction<3> *grad_d = new CompFunction<3>();
        for (int icomp = 0; icomp < inp.Ncomp(); icomp++) {
            for (int ocomp = 0; ocomp < 4; ocomp++) {
                if (std::norm(metric[icomp][ocomp]) > MachinePrec) {
                    grad_d->func_ptr->Ncomp = ocomp + 1;
                    if (inp.isreal()) {
                        grad_d->func_ptr->isreal = 1;
                        grad_d->func_ptr->iscomplex = 0;
                        grad_d->CompD[ocomp] = new FunctionTree<3, double>(inp.CompD[0]->getMRA());
                        apply(*(grad_d->CompD[ocomp]), oper, *inp.CompD[icomp], d);
                        if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                            grad_d->CompD[ocomp]->rescale((metric[icomp][ocomp]).real());
                        }
                    } else {
                        grad_d->func_ptr->isreal = 0;
                        grad_d->func_ptr->iscomplex = 1;
                        grad_d->CompC[ocomp] = new FunctionTree<3, ComplexDouble>(inp.CompC[0]->getMRA());
                        apply(*(grad_d->CompC[ocomp]), oper, *inp.CompC[icomp], d);
                        if (abs(metric[icomp][ocomp] - 1.0) > MachinePrec) {
                            grad_d->CompC[ocomp]->rescale(metric[icomp][ocomp]);
                        }
                    }
                }
            }
        }
        out.push_back(grad_d);
    }
    return out;
}

/**
 * @brief Compute the divergence of a vector field using a derivative operator.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[out] out Output scalar function.
 * @param[in]  oper Derivative operator applied to each component.
 * @param[in]  inp  Vector of D function components with coefficients.
 *
 * @details
 * Applies the derivative to each component along its matching direction and
 * sums the results on the **union grid** of the widened component grids.
 *
 * @note The length of `inp` must equal `D`. The output should contain only
 * empty root nodes on entry.
 */
template <int D, typename T>
void divergence(FunctionTree<D, T> &out,
                DerivativeOperator<D> &oper,
                FunctionTreeVector<D, T> &inp) {
    if (inp.size() != D) MSG_ABORT("Dimension mismatch");
    for (auto i = 0; i < inp.size(); i++)
        if (out.getMRA() != get_func(inp, i).getMRA()) MSG_ABORT("Incompatible MRA");

    FunctionTreeVector<D, T> tmp_vec;
    for (int d = 0; d < D; d++) {
        T coef_d = get_coef(inp, d);
        FunctionTree<D, T> &func_d = get_func(inp, d);
        auto *out_d = new FunctionTree<D, T>(func_d.getMRA());
        apply(*out_d, oper, func_d, d);
        tmp_vec.push_back(std::make_tuple(coef_d, out_d));
    }
    build_grid(out, tmp_vec);
    add(-1.0, out, tmp_vec, 0);
    clear(tmp_vec, true);
}

/**
 * @brief Divergence for multi-component inputs with metric (not implemented).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 */
template <int D, typename T>
void divergence(CompFunction<D> &out,
                DerivativeOperator<D> &oper,
                FunctionTreeVector<D, T> *inp,
                const ComplexDouble (*metric)[4]) {
    MSG_ABORT("not implemented");
}

/**
 * @brief Convenience overload: divergence from a list of unweighted component trees.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 */
template <int D, typename T>
void divergence(FunctionTree<D, T> &out,
                DerivativeOperator<D> &oper,
                std::vector<FunctionTree<D, T> *> &inp) {
    FunctionTreeVector<D, T> inp_vec;
    for (auto &t : inp) inp_vec.push_back({1.0, t});
    divergence(out, oper, inp_vec);
}

/* ---------- Explicit template instantiations ---------- */

template void apply<1, double>(double prec, FunctionTree<1, double> &out, ConvolutionOperator<1> &oper, FunctionTree<1, double> &inp, int maxIter, bool absPrec);
template void apply<2, double>(double prec, FunctionTree<2, double> &out, ConvolutionOperator<2> &oper, FunctionTree<2, double> &inp, int maxIter, bool absPrec);
template void apply<3, double>(double prec, FunctionTree<3, double> &out, ConvolutionOperator<3> &oper, FunctionTree<3, double> &inp, int maxIter, bool absPrec);
template void apply<1>(double prec, CompFunction<1> &out, ConvolutionOperator<1> &oper, const CompFunction<1> &inp, const ComplexDouble (*metric)[4], int maxIter = -1, bool absPrec = false);
template void apply<2>(double prec, CompFunction<2> &out, ConvolutionOperator<2> &oper, const CompFunction<2> &inp, const ComplexDouble (*metric)[4], int maxIter = -1, bool absPrec = false);
template void apply<3>(double prec, CompFunction<3> &out, ConvolutionOperator<3> &oper, const CompFunction<3> &inp, const ComplexDouble (*metric)[4], int maxIter = -1, bool absPrec = false);

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

template void apply<1, ComplexDouble>(double prec,
                                      FunctionTree<1, ComplexDouble> &out,
                                      ConvolutionOperator<1> &oper,
                                      FunctionTree<1, ComplexDouble> &inp,
                                      FunctionTreeVector<1, ComplexDouble> &precTrees,
                                      int maxIter,
                                      bool absPrec);
template void apply<2, ComplexDouble>(double prec,
                                      FunctionTree<2, ComplexDouble> &out,
                                      ConvolutionOperator<2> &oper,
                                      FunctionTree<2, ComplexDouble> &inp,
                                      FunctionTreeVector<2, ComplexDouble> &precTrees,
                                      int maxIter,
                                      bool absPrec);
template void apply<3, ComplexDouble>(double prec,
                                      FunctionTree<3, ComplexDouble> &out,
                                      ConvolutionOperator<3> &oper,
                                      FunctionTree<3, ComplexDouble> &inp,
                                      FunctionTreeVector<3, ComplexDouble> &precTrees,
                                      int maxIter,
                                      bool absPrec);

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

template void apply(CompFunction<3> &out, DerivativeOperator<3> &oper, CompFunction<3> &inp, int dir = -1, const ComplexDouble (*metric)[4]);

} // namespace mrcpp