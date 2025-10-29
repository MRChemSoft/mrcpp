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
 * @file grid.cpp
 * @brief Utilities for constructing, copying, clearing, and refining
 *        multiresolution grids and functions.
 *
 * @details
 * This module provides a unified set of routines for:
 *
 * - **Uniform grid construction** by splitting all leaves a fixed number of times.
 * - **Analytic-driven/adaptive grid construction** using a
 *   #mrcpp::RepresentableFunction as a splitter oracle.
 * - **Gaussian-expansion–driven grid construction** that places resolution
 *   according to Gaussian positions and exponents (supports periodic and
 *   non-periodic worlds).
 * - **Copying grids** (structure only) and **copying functions** (coefficients)
 *   between trees with the same #mrcpp::MultiResolutionAnalysis.
 * - **Clearing** coefficients on an existing grid without altering its topology.
 * - **Refining** an existing grid either uniformly, by precision-driven
 *   wavelet criteria, by another reference tree, or by an analytic function.
 *
 * All routines operate on #mrcpp::FunctionTree objects (and component-wise on
 * #mrcpp::CompFunction where relevant). Behind the scenes, they use
 * #mrcpp::TreeBuilder with different adaptors:
 *
 * - #mrcpp::SplitAdaptor: unconditional splitting.
 * - #mrcpp::WaveletAdaptor: split by wavelet-based precision criterion.
 * - #mrcpp::AnalyticAdaptor: split by analytic visibility/zero checks.
 * - #mrcpp::CopyAdaptor: split to match an existing tree structure.
 *
 * @note Unless otherwise stated, all "build_grid" functions **extend** the
 *       current grid of the output tree; they do not clear it first. Use
 *       #copy_grid when you want the output to match another grid exactly
 *       (it clears first).
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

/**
 * @brief Build an **empty** grid by uniform refinement.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[out] out    Output tree whose grid is refined.
 * @param[in]  scales Number of uniform refinement sweeps to apply.
 *
 * @details
 * Performs `scales` iterations of unconditional splitting on **all** current
 * leaf nodes (using #mrcpp::SplitAdaptor). No coefficients are created; this
 * only modifies the grid topology.
 *
 * @note Starts from the existing grid of @p out and extends it.
 */
template <int D, typename T> void build_grid(FunctionTree<D, T> &out, int scales) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    DefaultCalculator<D, T> calculator;
    SplitAdaptor<D, T> adaptor(maxScale, true); // Splits all nodes
    for (auto n = 0; n < scales; n++) builder.build(out, calculator, adaptor, 1);
}

/**
 * @brief Build an **empty** grid guided by an analytic function (adaptive).
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[out] out     Output tree whose grid will be extended.
 * @param[in]  inp     Analytic function used as a splitting oracle.
 * @param[in]  maxIter Maximum number of refinement iterations (-1 = unbounded).
 *
 * @details
 * Uses #mrcpp::AnalyticAdaptor to ask the analytic function @p inp whether a
 * node is visible at a given scale and whether it is identically zero on the
 * node interval. Nodes are split until convergence or @p maxIter is reached.
 *
 * @note Requires @p inp to implement `isVisibleAtScale()` and
 *       `isZeroOnInterval()`.
 */
template <int D, typename T> void build_grid(FunctionTree<D, T> &out, const RepresentableFunction<D, T> &inp, int maxIter) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    AnalyticAdaptor<D, T> adaptor(inp, maxScale);
    DefaultCalculator<D, T> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    print::separator(10, ' ');
}

/**
 * @brief Build an **empty** grid guided by a Gaussian expansion (adaptive).
 *
 * @tparam D Spatial dimension.
 * @param[out] out     Output tree whose grid will be extended.
 * @param[in]  inp     Gaussian expansion.
 * @param[in]  maxIter Maximum number of refinement iterations (-1 = unbounded).
 *
 * @details
 * For a non-periodic world:
 *   iterates over all Gaussians in @p inp and drives refinement with
 *   #mrcpp::AnalyticAdaptor using each Gaussian's position and exponent.
 *
 * For a periodic world:
 *   copies and reuses the same logic via temporary Gaussian objects so that
 *   periodic replication is handled consistently.
 *
 * Higher exponents imply finer resolution near the Gaussian center.
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
        auto period = out.getMRA().getWorldBox().getScalingFactors();
        (void)period; // currently unused; kept to document intent
        for (auto i = 0; i < inp.size(); i++) {
            auto *gauss = inp.getFunc(i).copy();
            build_grid(out, *gauss, maxIter);
            delete gauss;
        }
    }
    print::separator(10, ' ');
}

/**
 * @brief Build an **empty** grid by taking the union with another MW tree.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[out] out     Output tree to be extended.
 * @param[in]  inp     Input tree whose structure drives refinement.
 * @param[in]  maxIter Maximum number of refinement iterations (-1 = unbounded).
 *
 * @details
 * Uses #mrcpp::CopyAdaptor to ensure that any node that exists (and has
 * children) in @p inp will also exist in @p out after the call.
 *
 * @warning @p out and @p inp must share the same #mrcpp::MultiResolutionAnalysis.
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

/**
 * @brief Build an **empty** grid by taking the union of several MW trees.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[out] out     Output tree to be extended.
 * @param[in]  inp     Vector of (coef, tree) pairs.
 * @param[in]  maxIter Maximum number of refinement iterations (-1 = unbounded).
 *
 * @details
 * Uses #mrcpp::CopyAdaptor to extend @p out so that all nodes present in any
 * of the input trees are represented in the resulting grid (union).
 *
 * @warning All trees must share the same #mrcpp::MultiResolutionAnalysis as @p out.
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

/**
 * @brief Convenience overload: build a grid from a list of tree pointers.
 */
template <int D, typename T> void build_grid(FunctionTree<D, T> &out, std::vector<FunctionTree<D, T> *> &inp, int maxIter) {
    FunctionTreeVector<D, T> inp_vec;
    for (auto *t : inp) inp_vec.push_back({1.0, t});
    build_grid(out, inp_vec, maxIter);
}

/**
 * @brief Copy a function from one tree to the fixed grid of another.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[out] out Output tree (grid must already exist).
 * @param[in]  inp Input tree (source of coefficients).
 *
 * @details
 * Traverses the **current leaves** of @p out and copies the corresponding
 * coefficients from @p inp where nodes align, using the addition kernel
 * with fixed grid (no refinement).
 *
 * @note Overwrites existing coefficients in @p out; does not modify its grid.
 */
template <int D, typename T> void copy_func(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    FunctionTreeVector<D, T> tmp_vec;
    tmp_vec.push_back(std::make_tuple(1.0, &inp));
    add(-1.0, out, tmp_vec);
}

/**
 * @brief Make @p out's grid an exact copy of @p inp's grid (clears first).
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[out] out Output tree to be rebuilt.
 * @param[in]  inp Input tree supplying the grid structure.
 *
 * @details
 * Clears @p out completely (removes all nodes) and then extends its grid to
 * match @p inp using #build_grid(out, inp).
 *
 * @warning @p out and @p inp must share the same #mrcpp::MultiResolutionAnalysis.
 */
template <int D, typename T> void copy_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA")
    out.clear();
    build_grid(out, inp);
}

/**
 * @brief Component-wise grid copy for composite functions (clears first).
 *
 * @tparam D Spatial dimension.
 * @param[out] out Destination composite function.
 * @param[in]  inp Source composite function.
 *
 * @details
 * Recreates @p out with the same number of components and data parameters as
 * @p inp, then for each component copies the grid using the tree-based
 * #build_grid overload.
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

/**
 * @brief Clear coefficients on an existing grid (topology unchanged).
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[in,out] out Tree whose coefficients will be zeroed.
 *
 * @details
 * Uses #mrcpp::TreeBuilder::clear with #mrcpp::DefaultCalculator to reset
 * coefficients while preserving node structure.
 */
template <int D, typename T> void clear_grid(FunctionTree<D, T> &out) {
    TreeBuilder<D, T> builder;
    DefaultCalculator<D, T> calculator;
    builder.clear(out, calculator);
}

/**
 * @brief Uniformly refine a grid and **transfer scaling coefficients**.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[in,out] out    Tree to refine.
 * @param[in]     scales Number of refinement sweeps.
 * @return Number of nodes that were split.
 *
 * @details
 * Splits all leaves `scales` times using #mrcpp::TreeBuilder::split with
 * coefficient transfer to children, so the function representation remains
 * unchanged while resolution increases.
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

/**
 * @brief Precision-driven refinement using wavelet criteria.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[in,out] out     Tree to refine.
 * @param[in]     prec    Precision target for split checks.
 * @param[in]     absPrec If true, use absolute precision; otherwise relative.
 * @return Number of nodes that were split.
 *
 * @details
 * Uses #mrcpp::WaveletAdaptor to test split conditions based on wavelet
 * coefficients against @p prec (absolute or relative). When splitting, scales
 * are updated by transferring coefficients to the children.
 */
template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, double prec, bool absPrec) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    int nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

/**
 * @brief Refine a grid to include all structure present in a reference tree.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[in,out] out Tree to refine (and receive coefficient transfer).
 * @param[in]     inp Reference tree that defines where @p out should split.
 * @return Number of nodes that were split.
 *
 * @details
 * Uses #mrcpp::CopyAdaptor to mirror structural refinement from @p inp into
 * @p out and transfers coefficients to children where splits occur.
 */
template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, FunctionTree<D, T> &inp) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA")
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    CopyAdaptor<D, T> adaptor(inp, maxScale, nullptr);
    auto nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

/**
 * @brief Analytic-driven refinement using a representable function.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar coefficient type.
 * @param[in,out] out Tree to refine.
 * @param[in]     inp Analytic function to act as a split oracle.
 * @return Number of nodes that were split.
 *
 * @details
 * Uses #mrcpp::AnalyticAdaptor to request refinement where @p inp is visible
 * at scale and not identically zero on the cell. Coefficients are transferred
 * upon splitting so the represented function remains unchanged.
 */
template <int D, typename T> int refine_grid(FunctionTree<D, T> &out, const RepresentableFunction<D, T> &inp) {
    auto maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    AnalyticAdaptor<D, T> adaptor(inp, maxScale);
    int nSplit = builder.split(out, adaptor, true);
    return nSplit;
}

// -------------------- explicit instantiations --------------------

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
