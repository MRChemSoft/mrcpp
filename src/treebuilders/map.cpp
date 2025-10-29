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
 * @file map.cpp
 * @brief Adaptive mapping of multiresolution (MW) function trees through a user
 *        supplied scalar-to-scalar mapping.
 *
 * @details
 * This module implements an adaptive **pointwise mapping** of an input
 * #mrcpp::FunctionTree onto an output #mrcpp::FunctionTree by applying a user
 * provided mapping function \f$f:\mathbb{R}\to\mathbb{R}\f$ to the function
 * values represented on the MW grid.
 *
 * The mapping is realized via the standard MRCPP build loop:
 *  - On the **current** output grid, coefficients are computed by evaluating
 *    the input function and applying the mapping function (handled by
 *    #mrcpp::MapCalculator).
 *  - A **wavelet-based split criterion** (via #mrcpp::WaveletAdaptor) refines
 *    the grid wherever the mapped function requires more resolution to meet
 *    the requested precision.
 *  - This **refine–recompute** cycle repeats until convergence or a maximum
 *    number of iterations is reached.
 *
 * ### Precision semantics
 * - If `absPrec == false` (default), the adaptor uses **relative precision**:
 *   refinement stops when wavelet coefficients are small compared to the
 *   current function norm, roughly \f$|d| < \varepsilon\,/\,\|f\|\f$.
 * - If `absPrec == true`, the adaptor enforces an **absolute threshold**:
 *   \f$|d| < \varepsilon\f$.
 *
 * ### Responsibilities and caveats
 * - MRCPP does **not** impose constraints on the mapping function; the user
 *   must ensure it is numerically safe (no division by zero, no overflow, etc.).
 * - The mapping is **pointwise**: it does not solve PDEs or apply operators.
 *   For linear/nonlinear operators, consider specialized operator modules.
 *
 * ### Typical usage
 * @code
 * // Assume 'mra' is a configured MultiResolutionAnalysis<D>
 * FunctionTree<3,double> in(mra), out(mra);
 * // ... build 'in' somehow (project analytic function, read from file, etc.)
 *
 * auto clamp_nonnegative = [](double x) { return x < 0.0 ? 0.0 : x; };
 * map<3>(1e-6, out, in, clamp_nonnegative, -1 /* maxIter (unbounded) */, false /* relative precision */);
 * @endcode
 *
 * @see mrcpp::MapCalculator, mrcpp::WaveletAdaptor, mrcpp::TreeBuilder
 */

#include "map.h"
#include "MapCalculator.h"
#include "MultiplicationCalculator.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "grid.h"
#include "trees/FunctionNode.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include <Eigen/Core>

namespace mrcpp {

/**
 * @brief Adaptively map an input MW function through a scalar mapping function.
 *
 * @tparam D        Spatial dimension (1, 2, or 3).
 *
 * @param[in]  prec     Target build precision (relative or absolute depending on @p absPrec).
 * @param[out] out      Output function tree to be constructed (should start empty).
 * @param[in]  inp      Input function tree providing the source values.
 * @param[in]  fmap     Mapping function \f$f:\mathbb{R}\to\mathbb{R}\f$ to apply pointwise.
 * @param[in]  maxIter  Maximum refinement iterations (negative = unbounded).
 * @param[in]  absPrec  If true: interpret @p prec as absolute; otherwise relative.
 *
 * @details
 * Pipeline:
 *  1. Create a #mrcpp::MapCalculator that evaluates @p inp and applies @p fmap.
 *  2. Drive refinement with a #mrcpp::WaveletAdaptor at the MRA max scale,
 *     honoring @p prec and @p absPrec.
 *  3. Build the output via #mrcpp::TreeBuilder until convergence or @p maxIter.
 *  4. Perform bottom-up MW transform and square-norm computation for diagnostics.
 *  5. Clean temporary/generated artifacts on the input tree.
 *
 * @note
 * - The algorithm **extends** whatever grid @p out currently has. For a fresh build,
 *   ensure @p out is empty (no coefficients).
 * - The input and output trees must belong to a compatible MRA setup.
 *
 * @warning
 * The user is responsible for the numerical stability of @p fmap.
 * Discontinuous or extremely steep mappings may require tighter precision or
 * more iterations to resolve features adequately.
 */
template <int D>
void map(double prec,
         FunctionTree<D, double> &out,
         FunctionTree<D, double> &inp,
         FMap<double, double> fmap,
         int maxIter,
         bool absPrec) {

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, double> builder;
    WaveletAdaptor<D, double> adaptor(prec, maxScale, absPrec);
    MapCalculator<D, double> calculator(fmap, inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    inp.deleteGenerated();
    clean_t.stop();

    print::time(10, "Time transform", trans_t);
    print::time(10, "Time cleaning", clean_t);
    print::separator(10, ' ');
}

// explicit instantiations
template void map<1>(double prec, FunctionTree<1, double> &out, FunctionTree<1, double> &inp, FMap<double, double> fmap, int maxIter, bool absPrec);
template void map<2>(double prec, FunctionTree<2, double> &out, FunctionTree<2, double> &inp, FMap<double, double> fmap, int maxIter, bool absPrec);
template void map<3>(double prec, FunctionTree<3, double> &out, FunctionTree<3, double> &inp, FMap<double, double> fmap, int maxIter, bool absPrec);

} // Namespace mrcpp
