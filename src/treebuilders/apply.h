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

#pragma once
/**
 * @file apply.h
 * @brief Adaptive application of convolution/derivative operators to
 *        multiresolution (MW) functions and composite (multi-component) functions.
 *
 * @details
 * This header declares a family of routines that:
 * - apply **separable convolution operators** (near-/far-field or full) to MW trees,
 * - apply **derivative operators** to scalar or vector fields,
 * - compute **divergence** of vector fields, and
 * - compute **gradients**.
 *
 * Overloads exist for scalar MW trees (`FunctionTree<D,T>`) and for composite
 * multi-component fields (`CompFunction<D>`). For composite variants a
 * 4×4 complex **metric** can be supplied to define the componentwise inner
 * product / mixing; by default the identity metric is used.
 *
 * Precision and adaptivity:
 * - `prec` is the target build precision used by the adaptive refinement loop.
 * - `absPrec = false` → relative criterion; `true` → absolute threshold.
 * - `maxIter < 0` removes the iteration cap.
 */

#include "trees/FunctionTreeVector.h"
#include "utils/CompFunction.h"

namespace mrcpp {

// clang-format off

template <int D, typename T> class FunctionTree;
template <int D> class DerivativeOperator;
template <int D> class ConvolutionOperator;

/**
 * @brief Default 4×4 complex metric (identity).
 *
 * @details
 * Used by composite-function overloads to define component coupling / inner
 * product when applying operators. The default is the identity, i.e., no
 * cross-component mixing.
 */
constexpr ComplexDouble  defaultMetric [4][4] ={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

/**
 * @name Convolution application (scalar FunctionTree)
 * @{
 */

/**
 * @brief Apply a separable convolution operator adaptively.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (e.g., double, ComplexDouble).
 *
 * @param[in]  prec     Target precision for the adaptive build.
 * @param[out] out      Output function tree (built/extended adaptively).
 * @param[in]  oper     Convolution operator to apply.
 * @param[in]  inp      Input function tree.
 * @param[in]  maxIter  Maximum refinement iterations (-1 = unbounded).
 * @param[in]  absPrec  Use absolute (true) or relative (false) precision.
 */
template <int D, typename T>
void apply(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper,
           FunctionTree<D, T> &inp, int maxIter = -1, bool absPrec = false);

/**
 * @brief Apply a convolution operator with **per-node precision modulation**.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  prec       Base precision.
 * @param[out] out        Output function tree.
 * @param[in]  oper       Convolution operator.
 * @param[in]  inp        Input function tree.
 * @param[in]  precTrees  Vector of trees used to modulate local precision
 *                        (e.g., via node-wise scaling factors).
 * @param[in]  maxIter    Maximum refinement iterations.
 * @param[in]  absPrec    Absolute vs. relative precision.
 */
template <int D, typename T>
void apply(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper,
           FunctionTree<D, T> &inp, FunctionTreeVector<D, T> &precTrees,
           int maxIter = -1, bool absPrec = false);

/**
 * @brief Apply only the **far-field** contribution of a convolution operator.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output function tree.
 * @param[in]  oper     Convolution operator (far-field path will be used).
 * @param[in]  inp      Input function tree.
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute vs. relative precision.
 */
template <int D, typename T>
void apply_far_field(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper,
                     FunctionTree<D, T> &inp, int maxIter = -1, bool absPrec = false);

/**
 * @brief Apply only the **near-field** contribution of a convolution operator.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output function tree.
 * @param[in]  oper     Convolution operator (near-field path will be used).
 * @param[in]  inp      Input function tree.
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute vs. relative precision.
 */
template <int D, typename T>
void apply_near_field(double prec, FunctionTree<D, T> &out, ConvolutionOperator<D> &oper,
                      FunctionTree<D, T> &inp, int maxIter = -1, bool absPrec = false);

/** @} */

/**
 * @name Convolution application (composite CompFunction)
 * @{
 */

/**
 * @brief Apply a convolution operator to a composite function with a metric.
 *
 * @tparam D Spatial dimension.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output composite function.
 * @param[in]  oper     Convolution operator.
 * @param[in]  inp      Input composite function.
 * @param[in]  metric   Optional 4×4 complex metric (defaults to identity).
 * @param[in]  maxIter  Maximum refinement iterations (-1 = unbounded).
 * @param[in]  absPrec  Absolute vs. relative precision.
 *
 * @note Components can be coupled via @p metric during accumulation.
 */
template <int D>
void apply(double prec, CompFunction<D> &out, ConvolutionOperator<D> &oper,
           const CompFunction<D> &inp, const ComplexDouble (*metric)[4] = defaultMetric,
           int maxIter = -1, bool absPrec = false);

/**
 * @brief Apply a convolution operator to a composite function with
 *        precision-modulating trees and optional metric.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type used in @p precTrees.
 *
 * @param[in]  prec       Base precision.
 * @param[out] out        Output composite function.
 * @param[in]  oper       Convolution operator.
 * @param[in]  inp        Input composite function.
 * @param[in]  precTrees  Optional per-node precision modulators (may be nullptr).
 * @param[in]  metric     Optional 4×4 complex metric (may be nullptr → identity).
 * @param[in]  maxIter    Maximum refinement iterations.
 * @param[in]  absPrec    Absolute vs. relative precision.
 */
template <int D, typename T>
void apply(double prec, CompFunction<D> &out, ConvolutionOperator<D> &oper,
           CompFunction<D> &inp, FunctionTreeVector<D, T> *precTrees,
           ComplexDouble (*metric)[4] = nullptr, int maxIter = -1, bool absPrec = false);

/**
 * @brief Apply only the **far-field** part to a composite function.
 *
 * @tparam D Spatial dimension.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output composite function.
 * @param[in]  oper     Convolution operator.
 * @param[in]  inp      Input composite function.
 * @param[in]  metric   Optional 4×4 complex metric (defaults to identity).
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute vs. relative precision.
 */
template <int D>
void apply_far_field(double prec, CompFunction<D> &out, ConvolutionOperator<D> &oper,
                     CompFunction<D> &inp, const ComplexDouble (*metric)[4] = defaultMetric,
                     int maxIter = -1, bool absPrec = false);

/**
 * @brief Apply only the **near-field** part to a composite function.
 *
 * @tparam D Spatial dimension.
 *
 * @param[in]  prec     Target precision.
 * @param[out] out      Output composite function.
 * @param[in]  oper     Convolution operator.
 * @param[in]  inp      Input composite function.
 * @param[in]  metric   Optional 4×4 complex metric (defaults to identity).
 * @param[in]  maxIter  Maximum refinement iterations.
 * @param[in]  absPrec  Absolute vs. relative precision.
 */
template <int D>
void apply_near_field(double prec, CompFunction<D> &out, ConvolutionOperator<D> &oper,
                      CompFunction<D> &inp, const ComplexDouble (*metric)[4] = defaultMetric,
                      int maxIter = -1, bool absPrec = false);

/** @} */

/**
 * @name Derivative application
 * @{
 */

/**
 * @brief Apply a derivative operator to a scalar MW function.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[out] out  Output tree (derivative result).
 * @param[in]  oper Derivative operator.
 * @param[in]  inp  Input function.
 * @param[in]  dir  Application direction (0..D-1). If negative, use the
 *                  operator’s internal direction.
 */
template <int D, typename T>
void apply(FunctionTree<D, T> &out, DerivativeOperator<D> &oper,
           FunctionTree<D, T> &inp, int dir = -1);

/**
 * @brief Apply a derivative operator to a composite function with a metric.
 *
 * @tparam D Spatial dimension.
 *
 * @param[out] out    Output composite function.
 * @param[in]  oper   Derivative operator.
 * @param[in]  inp    Input composite function.
 * @param[in]  dir    Application direction (0..D-1). If negative, use operator’s default.
 * @param[in]  metric Optional 4×4 complex metric (defaults to identity).
 */
template <int D>
void apply(CompFunction<D> &out, DerivativeOperator<D> &oper,
           CompFunction<D> &inp, int dir = -1,
           const ComplexDouble (*metric)[4] = defaultMetric);

/** @} */

/**
 * @name Divergence
 * @{
 */

/**
 * @brief Divergence of a vector field given as separate component trees.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[out] out  Output scalar field (divergence).
 * @param[in]  oper Derivative operator (used per direction).
 * @param[in]  inp  Vector of component trees (size D expected).
 */
template <int D, typename T>
void divergence(FunctionTree<D, T> &out, DerivativeOperator<D> &oper,
                FunctionTreeVector<D, T> &inp);

/**
 * @brief Divergence of a composite vector field with metric.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type used by the composite.
 *
 * @param[out] out    Output scalar composite function (divergence).
 * @param[in]  oper   Derivative operator (used per direction).
 * @param[in]  inp    Pointer to vector of component composite functions.
 * @param[in]  metric Optional 4×4 complex metric.
 */
template <int D, typename T>
void divergence(CompFunction<D> &out, DerivativeOperator<D> &oper,
                FunctionTreeVector<D, T> *inp,
                const ComplexDouble (*metric)[4] = defaultMetric);

/**
 * @brief Divergence of a vector field given as a raw list of component pointers.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[out] out  Output scalar field (divergence).
 * @param[in]  oper Derivative operator.
 * @param[in]  inp  Vector of pointers to component trees (size D expected).
 */
template <int D, typename T>
void divergence(FunctionTree<D, T> &out, DerivativeOperator<D> &oper,
                std::vector<FunctionTree<D, T> *> &inp);

/**
 * @brief Divergence for composite fields given as raw component pointers with metric.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type used by components.
 *
 * @param[out] out    Output scalar composite function.
 * @param[in]  oper   Derivative operator.
 * @param[in]  inp    Pointer to vector of component tree pointers.
 * @param[in]  metric Optional 4×4 complex metric.
 */
template <int D, typename T>
void divergence(CompFunction<D> &out, DerivativeOperator<D> &oper,
                std::vector<FunctionTree<D, T> *> *inp,
                const ComplexDouble (*metric)[4]  = defaultMetric);

/** @} */

/**
 * @name Gradient
 * @{
 */

/**
 * @brief Gradient of a scalar field (returns D component trees).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 *
 * @param[in]  oper Derivative operator (used per direction).
 * @param[in]  inp  Input scalar field.
 * @return Vector of D component trees with directional derivatives.
 */
template <int D, typename T>
FunctionTreeVector<D, T> gradient(DerivativeOperator<D> &oper, FunctionTree<D, T> &inp);

/**
 * @brief Gradient (3D) for composite fields, returning heap-allocated components.
 *
 * @param[in]  oper   3D derivative operator.
 * @param[in]  inp    Input composite function.
 * @param[in]  metric Optional 4×4 complex metric (defaults to identity).
 * @return Vector of pointers to newly allocated component composite functions
 *         representing the gradient. The caller owns and must delete them.
 */
std::vector<CompFunction<3>*> gradient(DerivativeOperator<3> &oper, CompFunction<3> &inp,
                                       const ComplexDouble  (*metric)[4] = defaultMetric);
// clang-format on

} // namespace mrcpp