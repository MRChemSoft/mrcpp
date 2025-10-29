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
 * @file add.cpp
 * @brief Adaptive summation of multiwavelet (MW) function trees.
 *
 * @details
 * This module provides a family of `add` routines that assemble the linear
 * combination of one or more MW functions into an output MW function on an
 * adaptively refined grid.
 *
 * The summation is performed by the generic @ref TreeBuilder orchestrating:
 * - an @ref AdditionCalculator that evaluates the local sum of input trees
 *   with their numerical coefficients (and optional complex conjugation),
 * - a @ref WaveletAdaptor that refines the output grid where needed to meet
 *   the requested precision.
 *
 * The core algorithm (all overloads):
 * - Compute MW coefficients of the sum on the **current** output grid.
 * - Refine the grid according to the precision target.
 * - Repeat until convergence or until a maximum number of refinement
 *   iterations is reached.
 * - Finally transform the output to the MW domain and compute its squared norm.
 *
 * Precision and iteration controls:
 * - `prec < 0` or `maxIter = 0` disables refinement (single pass on
 *   the existing output grid).
 * - `maxIter < 0` removes the iteration limit and refines until the
 *   precision criterion is satisfied.
 * - `absPrec = true` interprets `prec` as an absolute tolerance, otherwise
 *   it is treated as a relative criterion.
 *
 * Requirements:
 * - All input trees must share the same @ref MultiResolutionAnalysis as the
 *   output tree, otherwise the routine aborts.
 *
 * Notes:
 * - The routine starts from whatever grid is already present in `out`. This
 *   grid is expected to be empty in terms of coefficients.
 * - Generated nodes present in input trees are removed at the end (cleanup).
 */

#include <tuple>
#include <vector>

#include "AdditionCalculator.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/**
 * @brief Sum two MW functions (with scalar weights) into an output tree using adaptive refinement.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (e.g., double or ComplexDouble).
 *
 * @param[in]  prec      Target build precision for the output function.
 * @param[out] out       Output function tree to build (its current grid is used as a starting point).
 * @param[in]  a         Numerical coefficient multiplying `inp_a`.
 * @param[in]  inp_a     First input function tree.
 * @param[in]  b         Numerical coefficient multiplying `inp_b`.
 * @param[in]  inp_b     Second input function tree.
 * @param[in]  maxIter   Maximum number of refinement iterations.
 *                       Use a negative value to allow unbounded refinement.
 *                       Use zero to disable refinement (single-pass build).
 * @param[in]  absPrec   If true, interpret `prec` as an absolute tolerance;
 *                       otherwise interpret it as relative.
 * @param[in]  conjugate When `T` is complex, conjugate all input trees before summation.
 *
 * @details
 * Builds `out ≈ a * inp_a (+) b * inp_b` to the requested precision on an adaptively
 * refined grid. After the build, `out` is transformed to the MW domain and its squared
 * norm is computed. The input trees are not modified except that any generated nodes
 * created temporarily during the build are cleaned up.
 */
template <int D, typename T>
void add(double prec,
         FunctionTree<D, T> &out,
         T a, FunctionTree<D, T> &inp_a,
         T b, FunctionTree<D, T> &inp_b,
         int maxIter,
         bool absPrec,
         bool conjugate) {
    FunctionTreeVector<D, T> tmp_vec;
    tmp_vec.push_back(std::make_tuple(a, &inp_a));
    tmp_vec.push_back(std::make_tuple(b, &inp_b));
    add(prec, out, tmp_vec, maxIter, absPrec, conjugate);
}

/**
 * @brief Sum a vector of MW functions (with scalar weights) into an output tree using adaptive refinement.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (e.g., double or ComplexDouble).
 *
 * @param[in]  prec      Target build precision for the output function.
 * @param[out] out       Output function tree to build (its current grid is used as a starting point).
 * @param[in]  inp       Vector of pairs (weight, pointer-to-tree) to be summed.
 * @param[in]  maxIter   Maximum number of refinement iterations.
 *                       Use a negative value to allow unbounded refinement.
 *                       Use zero to disable refinement (single-pass build).
 * @param[in]  absPrec   If true, interpret `prec` as an absolute tolerance;
 *                       otherwise interpret it as relative.
 * @param[in]  conjugate When `T` is complex, conjugate all input trees before summation.
 *
 * @details
 * Builds `out ≈ Σ_i w_i * f_i` to the requested precision on an adaptively refined grid.
 * The routine:
 * - verifies that all inputs share the same MRA as `out`,
 * - constructs a @ref WaveletAdaptor with the precision policy,
 * - uses an @ref AdditionCalculator to evaluate the local sums,
 * - runs @ref TreeBuilder to refine and assemble,
 * - finishes with MW transform and squared norm computation,
 * - and finally deletes any generated nodes from inputs.
 */
template <int D, typename T>
void add(double prec,
         FunctionTree<D, T> &out,
         FunctionTreeVector<D, T> &inp,
         int maxIter,
         bool absPrec,
         bool conjugate) {
    for (auto i = 0; i < inp.size(); i++)
        if (out.getMRA() != get_func(inp, i).getMRA()) MSG_ABORT("Incompatible MRA");

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    AdditionCalculator<D, T> calculator(inp, conjugate);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();

    trans_t.stop();

    Timer clean_t;
    for (int i = 0; i < inp.size(); i++) {
        FunctionTree<D, T> &tree = get_func(inp, i);
        tree.deleteGenerated();
    }
    clean_t.stop();

    print::time(10, "Time transform", trans_t);
    print::time(10, "Time cleaning", clean_t);
    print::separator(10, ' ');
}

/**
 * @brief Convenience overload: sum a list of unweighted trees (weights set to 1).
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (e.g., double or ComplexDouble).
 *
 * @param[in]  prec      Target build precision for the output function.
 * @param[out] out       Output function tree.
 * @param[in]  inp       Vector of pointers to input trees (all weights taken as 1).
 * @param[in]  maxIter   Maximum number of refinement iterations (see other overload).
 * @param[in]  absPrec   Absolute-vs-relative precision flag.
 * @param[in]  conjugate Conjugate complex inputs before summation.
 *
 * @details
 * Internally wraps the list into a @ref FunctionTreeVector with unit weights
 * and forwards to the vector-based overload.
 */
template <int D, typename T>
void add(double prec,
         FunctionTree<D, T> &out,
         std::vector<FunctionTree<D, T> *> &inp,
         int maxIter,
         bool absPrec,
         bool conjugate) {
    FunctionTreeVector<D, T> inp_vec;
    for (auto &t : inp) inp_vec.push_back({1.0, t});
    add(prec, out, inp_vec, maxIter, absPrec, conjugate);
}

/* ------- Explicit template instantiations (double) ------- */
template void
add<1, double>(double prec, FunctionTree<1, double> &out, double a, FunctionTree<1, double> &tree_a, double b, FunctionTree<1, double> &tree_b, int maxIter, bool absPrec, bool conjugate);
template void
add<2, double>(double prec, FunctionTree<2, double> &out, double a, FunctionTree<2, double> &tree_a, double b, FunctionTree<2, double> &tree_b, int maxIter, bool absPrec, bool conjugate);
template void
add<3, double>(double prec, FunctionTree<3, double> &out, double a, FunctionTree<3, double> &tree_a, double b, FunctionTree<3, double> &tree_b, int maxIter, bool absPrec, bool conjugate);

template void add<1, double>(double prec, FunctionTree<1, double> &out, FunctionTreeVector<1, double> &inp, int maxIter, bool absPrec, bool conjugate);
template void add<2, double>(double prec, FunctionTree<2, double> &out, FunctionTreeVector<2, double> &inp, int maxIter, bool absPrec, bool conjugate);
template void add<3, double>(double prec, FunctionTree<3, double> &out, FunctionTreeVector<3, double> &inp, int maxIter, bool absPrec, bool conjugate);

template void add<1, double>(double prec, FunctionTree<1, double> &out, std::vector<FunctionTree<1, double> *> &inp, int maxIter, bool absPrec, bool conjugate);
template void add<2, double>(double prec, FunctionTree<2, double> &out, std::vector<FunctionTree<2, double> *> &inp, int maxIter, bool absPrec, bool conjugate);
template void add<3, double>(double prec, FunctionTree<3, double> &out, std::vector<FunctionTree<3, double> *> &inp, int maxIter, bool absPrec, bool conjugate);

/* ------- Explicit template instantiations (ComplexDouble) ------- */
template void add<1, ComplexDouble>(double prec,
                                    FunctionTree<1, ComplexDouble> &out,
                                    ComplexDouble a,
                                    FunctionTree<1, ComplexDouble> &tree_a,
                                    ComplexDouble b,
                                    FunctionTree<1, ComplexDouble> &tree_b,
                                    int maxIter,
                                    bool absPrec,
                                    bool conjugate);
template void add<2, ComplexDouble>(double prec,
                                    FunctionTree<2, ComplexDouble> &out,
                                    ComplexDouble a,
                                    FunctionTree<2, ComplexDouble> &tree_a,
                                    ComplexDouble b,
                                    FunctionTree<2, ComplexDouble> &tree_b,
                                    int maxIter,
                                    bool absPrec,
                                    bool conjugate);
template void add<3, ComplexDouble>(double prec,
                                    FunctionTree<3, ComplexDouble> &out,
                                    ComplexDouble a,
                                    FunctionTree<3, ComplexDouble> &tree_a,
                                    ComplexDouble b,
                                    FunctionTree<3, ComplexDouble> &tree_b,
                                    int maxIter,
                                    bool absPrec,
                                    bool conjugate);

template void add<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, FunctionTreeVector<1, ComplexDouble> &inp, int maxIter, bool absPrec, bool conjugate);
template void add<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, FunctionTreeVector<2, ComplexDouble> &inp, int maxIter, bool absPrec, bool conjugate);
template void add<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, FunctionTreeVector<3, ComplexDouble> &inp, int maxIter, bool absPrec, bool conjugate);

template void add<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, std::vector<FunctionTree<1, ComplexDouble> *> &inp, int maxIter, bool absPrec, bool conjugate);
template void add<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, std::vector<FunctionTree<2, ComplexDouble> *> &inp, int maxIter, bool absPrec, bool conjugate);
template void add<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, std::vector<FunctionTree<3, ComplexDouble> *> &inp, int maxIter, bool absPrec, bool conjugate);

} // namespace mrcpp