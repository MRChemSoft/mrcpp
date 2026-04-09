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
 * @file
 * @brief Grid construction, copying, clearing, and refinement helpers for multiresolution trees.
 *
 * @details
 * This header declares a family of utilities to *construct* and *modify* the
 * topology (grid) of @ref mrcpp::FunctionTree without necessarily computing or
 * moving coefficients. The functions support several sources:
 * analytic/representable functions, existing trees, vectors of trees, and
 * explicit scale counts.
 *
 * ### Conventions
 * - `D` is the spatial dimension (typically 1–3).
 * - `T` is the coefficient scalar type (`double` or `std::complex<double>`).
 * - Functions named `build_grid` create (or enlarge) the *tree structure*
 *   of `out` to be adequate for representing the given input(s).
 * - Functions named `copy_grid` copy only the *structure* (no coefficients).
 * - `copy_func` copies *both* structure and coefficients.
 * - Functions named `refine_grid` add resolution either explicitly by a scale
 *   count or adaptively by a precision criterion.
 * - Functions return `int` indicate the number of newly created end-nodes
 *   (i.e., how many refinements were actually performed).
 */

#include "functions/RepresentableFunction.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"
#include "utils/CompFunction.h"

namespace mrcpp {

/**
 * @brief Create a uniform grid of fixed depth.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[out] out Target tree whose topology will be (re)built.
 * @param[in]  scales Number of refinement steps from the current state.
 *
 * @details
 * Starting from the current `out` topology (typically roots), subdivide each
 * active end-node `scales` times so that a regular grid of depth increased by
 * `scales` is obtained. No coefficients are computed or modified.
 */
template <int D, typename T>
void build_grid(FunctionTree<D, T> &out, int scales);

/**
 * @brief Build an adaptive grid suitable for a Gaussian expansion.
 *
 * @tparam D Spatial dimension.
 * @param[out] out Target **real** tree to receive the grid.
 * @param[in]  inp Analytic Gaussian expansion used as refinement oracle.
 * @param[in]  maxIter Maximum refinement passes; negative means “unbounded”
 *                     until convergence by the internal criterion.
 *
 * @details
 * Iteratively refines the tree so that the structure can represent `inp`
 * within the library’s default per-node criterion (e.g., band-limited model or
 * local projection error). Coefficients are not guaranteed to be written.
 */
template <int D>
void build_grid(FunctionTree<D> &out, const GaussExp<D> &inp, int maxIter = -1);

/**
 * @brief Build an adaptive grid for a generic representable function.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[out] out Target tree.
 * @param[in]  inp Representable function serving as refinement oracle.
 * @param[in]  maxIter Maximum refinement passes; negative means unbounded.
 *
 * @details
 * Uses evaluations/projections of `inp` to determine where refinement is
 * needed so that the resulting grid can capture `inp` with the library’s
 * default tolerance heuristic.
 */
template <int D, typename T>
void build_grid(FunctionTree<D, T> &out, const RepresentableFunction<D, T> &inp, int maxIter = -1);

/**
 * @brief Build a grid that can accommodate another tree’s resolution/support.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[out] out Target tree to be enlarged/refined.
 * @param[in]  inp Source tree whose structure (support + finest scales) is used.
 *
 * @details
 * Ensures that `out` has at least the resolution present in `inp` wherever
 * `inp` has support (a *grid union* operation). Coefficients are not copied.
 */
template <int D, typename T>
void build_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp, int maxIter = -1);

/**
 * @brief Build a grid that is a union of a vector of trees.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[out] out Target tree.
 * @param[in]  inp Vector of trees whose supports/resolutions are merged.
 * @param[in]  maxIter Optional iteration cap for staged refinement strategies.
 */
template <int D, typename T>
void build_grid(FunctionTree<D, T> &out, FunctionTreeVector<D, T> &inp, int maxIter = -1);

/**
 * @brief Build a grid that is a union of a list of tree pointers.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[out] out Target tree.
 * @param[in]  inp List of tree pointers to merge.
 * @param[in]  maxIter Optional iteration cap for staged refinement strategies.
 */
template <int D, typename T>
void build_grid(FunctionTree<D, T> &out, std::vector<FunctionTree<D, T> *> &inp, int maxIter = -1);

/**
 * @brief Deep copy a tree structure *and* coefficients.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[out] out Destination tree (reallocated as needed).
 * @param[in]  inp Source tree.
 */
template <int D, typename T>
void copy_func(FunctionTree<D, T> &out, FunctionTree<D, T> &inp);

/**
 * @brief Copy only the tree topology (no coefficients).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[out] out Destination tree (structure rebuilt).
 * @param[in]  inp Source tree whose topology is replicated.
 */
template <int D, typename T>
void copy_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp);

/**
 * @brief Copy only the topology for all components of a composite function.
 *
 * @tparam D Spatial dimension.
 * @param[out] out Destination composite function (components allocated as needed).
 * @param[in]  inp Source composite function.
 *
 * @details
 * For each component present in @p inp, ensure @p out has a corresponding
 * component with identical tree structure. Coefficients are not copied.
 */
template <int D>
void copy_grid(CompFunction<D> &out, CompFunction<D> &inp);

/**
 * @brief Clear the grid topology (prune to roots, drop nodes).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[out] out Tree to clear; MRA association remains intact.
 */
template <int D, typename T>
void clear_grid(FunctionTree<D, T> &out);

/**
 * @brief Refine uniformly by a fixed number of scales.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[in,out] out Tree to refine.
 * @param[in]     scales Number of subdivision steps to apply.
 * @return Number of new end-nodes created by the refinement.
 */
template <int D, typename T>
int refine_grid(FunctionTree<D, T> &out, int scales);

/**
 * @brief Adaptive refinement driven by a precision threshold.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[in,out] out Tree to refine.
 * @param[in]     prec Target precision (threshold).
 * @param[in]     absPrec If `true`, interpret @p prec as absolute tolerance;
 *                        otherwise relative to a norm estimate.
 * @return Number of new end-nodes created.
 *
 * @details
 * Subdivides those nodes whose local error/indicator exceeds the requested
 * threshold. The precise indicator depends on the library configuration
 * (e.g., wavelet-norm-based splitting).
 */
template <int D, typename T>
int refine_grid(FunctionTree<D, T> &out, double prec, bool absPrec = false);

/**
 * @brief Refine `out` so that its grid is at least as fine as `inp`.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[in,out] out Destination tree to refine.
 * @param[in]     inp Source tree providing the target finest scales.
 * @return Number of new end-nodes created.
 */
template <int D, typename T>
int refine_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp);

/**
 * @brief Adaptive refinement using a representable function as oracle.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param[in,out] out Tree to refine.
 * @param[in]     inp Representable function guiding refinement.
 * @return Number of new end-nodes created.
 *
 * @details
 * Samples or projects @p inp on candidate nodes and refines where the
 * estimated local error is above the internal criterion, creating a grid
 * appropriate for subsequently projecting @p inp.
 */
template <int D, typename T>
int refine_grid(FunctionTree<D, T> &out, const RepresentableFunction<D, T> &inp);

} // namespace mrcpp