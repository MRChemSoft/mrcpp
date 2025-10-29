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
 * @file project.cpp
 * @brief Projection of analytic (scalar or vector) functions onto a
 *        multiwavelet (MW) basis on an adaptively refined grid.
 *
 * @details
 * This module builds a MW representation of an analytic function by
 * adaptively refining the grid and computing (scale-/wavelet-) coefficients
 * until a user-prescribed tolerance is achieved.
 *
 * ### Algorithm (adaptive projection)
 * 1. Start from the current grid in @p out (should be empty or root-only).
 * 2. On the current leaves, compute MW coefficients using
 *    ProjectionCalculator (quadrature in the scaling basis).
 * 3. Use WaveletAdaptor to decide where to refine:
 *    - **Relative precision** (default): stop when local wavelet norms
 *      drop below `prec * ||f||_node`.
 *    - **Absolute precision** (`absPrec = true`): stop when local wavelet
 *      norms drop below `prec`.
 * 4. Repeat until convergence or `maxIter` is reached.
 * 5. Perform final MW transforms (TopDown/BottomUp as needed) and compute
 *    the tree square-norm for bookkeeping.
 *
 * The projection accounts for non-unit world-box scaling through
 * a per-dimension scaling factor passed to ProjectionCalculator.
 *
 * @note The functions here operate on templated dimension @p D (1,2,3)
 *       and coefficient type @p T (double or ComplexDouble).
 */

#include "project.h"
#include "ProjectionCalculator.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "functions/AnalyticFunction.h"
#include "trees/FunctionTree.h"
#include "trees/MultiResolutionAnalysis.h"

#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/**
 * @brief Project a lambda/std::function onto the MW basis (convenience overload).
 *
 * Wraps the callable into an AnalyticFunction and delegates to the
 * RepresentableFunction overload.
 *
 * @tparam D Spatial dimension (1,2,3).
 * @tparam T Coefficient type (double or ComplexDouble).
 *
 * @param[in] prec    Target precision (relative by default, see @p absPrec).
 * @param[out] out    Output function tree to be built (should contain only empty roots).
 * @param[in] func    Callable \f$f:\mathbb{R}^D \to T\f$ returning values at coordinates.
 * @param[in] maxIter Maximum refinement iterations (-1 = no bound).
 * @param[in] absPrec Use absolute (true) or relative (false, default) thresholding.
 *
 * @details
 * This is syntactic sugar for quickly projecting a user-provided callable.
 * The adaptive procedure, grid policy, and stopping criteria are identical
 * to the main projection overload below.
 *
 * @note The current grid in @p out is honored and extended; it is not cleared.
 *       For a fresh build, ensure @p out has only root nodes and no coefficients.
 */
template <int D, typename T>
void project(double prec,
             FunctionTree<D, T> &out,
             std::function<T(const Coord<D> &r)> func,
             int maxIter,
             bool absPrec) {
    AnalyticFunction<D, T> inp(func);
    mrcpp::project(prec, out, inp, maxIter, absPrec);
}

/**
 * @brief Project a RepresentableFunction onto the MW basis, adaptive grid.
 *
 * @tparam D Spatial dimension (1,2,3).
 * @tparam T Coefficient type (double or ComplexDouble).
 *
 * @param[in] prec    Target precision (relative by default, see @p absPrec).
 * @param[out] out    Output function tree to be built (should contain only empty roots).
 * @param[in] inp     Analytic/representable function to project.
 * @param[in] maxIter Maximum number of refinement iterations (-1 = unbounded).
 * @param[in] absPrec Use absolute (true) or relative (false) thresholding.
 *
 * @details
 * - Builds a WaveletAdaptor with precision policy (relative/absolute).
 * - Creates a ProjectionCalculator configured with world-box scaling
 *   factors to ensure correct physical rescaling of integrals.
 * - Uses TreeBuilder to iterate:
 *      compute coefs → test refinement → split nodes → repeat.
 * - Finalizes with a BottomUp MW transform and tree norm accumulation.
 *
 * @par Precision semantics
 * - **Relative** (`absPrec=false`): local wavelet norm compared to local function norm.
 * - **Absolute** (`absPrec=true`): local wavelet norm compared to @p prec directly.
 *
 * @warning The output tree @p out must be compatible (same MRA/world box)
 *          with any other trees you later combine it with.
 */
template <int D, typename T>
void project(double prec,
             FunctionTree<D, T> &out,
             RepresentableFunction<D, T> &inp,
             int maxIter,
             bool absPrec) {
    int maxScale = out.getMRA().getMaxScale();
    const auto scaling_factor = out.getMRA().getWorldBox().getScalingFactors();

    TreeBuilder<D, T> builder;
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    ProjectionCalculator<D, T> calculator(inp, scaling_factor);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');
}

/**
 * @brief Project a vector of analytic functions (component-wise), adaptive grid.
 *
 * @tparam D Spatial dimension (1,2,3).
 * @tparam T Coefficient type (double or ComplexDouble).
 *
 * @param[in] prec    Target precision (relative by default, see @p absPrec).
 * @param[out] out    Output vector of trees (size must match @p func).
 * @param[in] func    Vector of component callables \f$f_j:\mathbb{R}^D \to T\f$.
 * @param[in] maxIter Maximum refinement iterations (-1 = unbounded).
 * @param[in] absPrec Use absolute (true) or relative (false) thresholding.
 *
 * @details
 * Projects each component independently with the same precision policy and
 * refinement limits, storing the result in the corresponding entry of @p out.
 *
 * @throws MSG_ABORT if @p out.size() != @p func.size().
 */
template <int D, typename T>
void project(double prec,
             FunctionTreeVector<D, T> &out,
             std::vector<std::function<T(const Coord<D> &r)>> func,
             int maxIter,
             bool absPrec) {
    if (out.size() != func.size()) MSG_ABORT("Size mismatch");
    for (auto j = 0; j < D; j++) {
        mrcpp::project<D>(prec, get_func(out, j), func[j], maxIter, absPrec);
    }
}

template void project<1, double>(double prec, FunctionTree<1, double> &out, RepresentableFunction<1, double> &inp, int maxIter, bool absPrec);
template void project<2, double>(double prec, FunctionTree<2, double> &out, RepresentableFunction<2, double> &inp, int maxIter, bool absPrec);
template void project<3, double>(double prec, FunctionTree<3, double> &out, RepresentableFunction<3, double> &inp, int maxIter, bool absPrec);

template void project<1, double>(double prec, FunctionTree<1, double> &out, std::function<double(const Coord<1> &r)> func, int maxIter, bool absPrec);
template void project<2, double>(double prec, FunctionTree<2, double> &out, std::function<double(const Coord<2> &r)> func, int maxIter, bool absPrec);
template void project<3, double>(double prec, FunctionTree<3, double> &out, std::function<double(const Coord<3> &r)> func, int maxIter, bool absPrec);
template void project<1, double>(double prec, FunctionTreeVector<1, double> &out, std::vector<std::function<double(const Coord<1> &r)>> inp, int maxIter, bool absPrec);
template void project<2, double>(double prec, FunctionTreeVector<2, double> &out, std::vector<std::function<double(const Coord<2> &r)>> inp, int maxIter, bool absPrec);
template void project<3, double>(double prec, FunctionTreeVector<3, double> &out, std::vector<std::function<double(const Coord<3> &r)>> inp, int maxIter, bool absPrec);

template void project<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, RepresentableFunction<1, ComplexDouble> &inp, int maxIter, bool absPrec);
template void project<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, RepresentableFunction<2, ComplexDouble> &inp, int maxIter, bool absPrec);
template void project<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, RepresentableFunction<3, ComplexDouble> &inp, int maxIter, bool absPrec);

template void project<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, std::function<ComplexDouble(const Coord<1> &r)> func, int maxIter, bool absPrec);
template void project<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, std::function<ComplexDouble(const Coord<2> &r)> func, int maxIter, bool absPrec);
template void project<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, std::function<ComplexDouble(const Coord<3> &r)> func, int maxIter, bool absPrec);
template void project<1, ComplexDouble>(double prec, FunctionTreeVector<1, ComplexDouble> &out, std::vector<std::function<ComplexDouble(const Coord<1> &r)>> inp, int maxIter, bool absPrec);
template void project<2, ComplexDouble>(double prec, FunctionTreeVector<2, ComplexDouble> &out, std::vector<std::function<ComplexDouble(const Coord<2> &r)>> inp, int maxIter, bool absPrec);
template void project<3, ComplexDouble>(double prec, FunctionTreeVector<3, ComplexDouble> &out, std::vector<std::function<ComplexDouble(const Coord<3> &r)>> inp, int maxIter, bool absPrec);

} // namespace mrcpp
