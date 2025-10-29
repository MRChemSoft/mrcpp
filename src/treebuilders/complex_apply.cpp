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
 * @file complex_apply.cpp
 * @brief Complex-valued application of multiresolution convolution operators.
 *
 * @details
 * This module provides a **complex** front-end to the real-valued adaptive
 * application pipeline used throughout MRCPP. A complex operator
 * \f$ \mathcal{O} = \mathcal{O}_\mathrm{R} + i\,\mathcal{O}_\mathrm{I} \f$
 * acting on a complex function
 * \f$ f = f_\mathrm{R} + i\,f_\mathrm{I} \f$
 * is evaluated via the standard decomposition:
 * \f[
 *   \mathcal{O} f
 *   = (\mathcal{O}_\mathrm{R} f_\mathrm{R} - \mathcal{O}_\mathrm{I} f_\mathrm{I})
 *     \;+\;
 *     i\,(\mathcal{O}_\mathrm{I} f_\mathrm{R} + \mathcal{O}_\mathrm{R} f_\mathrm{I}).
 * \f]
 *
 * Internally, the routine delegates every real application
 * \f$ \mathcal{O}_\bullet f_\bullet \f$ to the standard adaptive `apply` for
 * real data structures (see `apply.h`), and then combines the four real
 * results to produce the complex output.
 *
 * ### Precision model and adaptivity
 * The same adaptive refinement loop is honored as in the real case:
 * - **Relative precision** (default): refine where local wavelet details exceed
 *   a fraction of the local norm.
 * - **Absolute precision** (`absPrec = true`): refine until local details fall
 *   below a fixed absolute threshold.
 *
 * The `prec` parameter and `maxIter` semantics are identical to the real-valued
 * `apply`:
 * - `prec < 0` or `maxIter = 0` disables refinement,
 * - `maxIter < 0` removes the iteration bound.
 *
 * ### Preconditions
 * - All real and imaginary parts (operator and function) must share the same
 *   `MultiResolutionAnalysis`.
 * - The output complex object should reference **empty** (uninitialized) trees
 *   at entry; the routine will construct their contents.
 *
 * @note This is a thin complex wrapper; all heavy lifting (bandwidth computation,
 *       adaptive splitting, transformations, norm updates) happens in the
 *       underlying real `apply`.
 */

#include "complex_apply.h"
#include "ConvolutionCalculator.h"
#include "CopyAdaptor.h"
#include "DefaultCalculator.h"
#include "DerivativeCalculator.h"
#include "SplitAdaptor.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "apply.h"
#include "grid.h"
#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/**
 * @brief Apply a complex convolution operator to a complex function (adaptive).
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 *
 * @param[in]  prec     Target build precision for the adaptive application.
 * @param[out] out      Complex output function tree (real and imaginary parts filled).
 * @param[in]  oper     Complex convolution operator (real and imaginary parts provided).
 * @param[in]  inp      Complex input function tree (real and imaginary parts provided).
 * @param[in]  maxIter  Maximum refinement iterations; `-1` means unbounded.
 * @param[in]  absPrec  Use absolute (`true`) versus relative (`false`, default) precision.
 *
 * @details
 * The routine evaluates
 * \f[
 *   \Re(\mathcal{O}f) = \mathcal{O}_\mathrm{R} f_\mathrm{R} - \mathcal{O}_\mathrm{I} f_\mathrm{I},\quad
 *   \Im(\mathcal{O}f) = \mathcal{O}_\mathrm{I} f_\mathrm{R} + \mathcal{O}_\mathrm{R} f_\mathrm{I}
 * \f]
 * by two real `apply` calls per part, followed by linear combinations via `add`.
 * Temporary real trees are allocated on the same MRA as the input.
 *
 * ### Implementation notes
 * - The real building blocks `apply(prec, ...)` are identical to the scalar path
 *   and include: bandwidth precomputation, adaptive refinement, top-down coarse
 *   contributions, bottom-up transforms, and norm updates.
 * - Output parts are formed with `add(prec, ...)` to maintain consistent grid
 *   and transformation state.
 *
 * @warning The MRA of `inp.real`, `inp.imaginary`, `oper.real`, and
 *          `oper.imaginary` must match. No cross-MRA application is supported.
 */
template <int D>
void apply(double prec,
           ComplexObject<FunctionTree<D>> &out,
           ComplexObject<ConvolutionOperator<D>> &oper,
           ComplexObject<FunctionTree<D>> &inp,
           int maxIter,
           bool absPrec) {
    FunctionTree<D> temp1(inp.real->getMRA());
    FunctionTree<D> temp2(inp.real->getMRA());

    // Real part:  OR*FR - OI*FI
    apply(prec, temp1, *oper.real,      *inp.real,      maxIter, absPrec);
    apply(prec, temp2, *oper.imaginary, *inp.imaginary, maxIter, absPrec);
    add(prec, *out.real, 1.0, temp1, -1.0, temp2);

    // Imag part:  OI*FR + OR*FI
    apply(prec, temp1, *oper.imaginary, *inp.real,      maxIter, absPrec);
    apply(prec, temp2, *oper.real,      *inp.imaginary, maxIter, absPrec);
    add(prec, *out.imaginary, 1.0, temp1, 1.0, temp2);
}

template void apply<1>(double prec,
                       ComplexObject<FunctionTree<1>> &out,
                       ComplexObject<ConvolutionOperator<1>> &oper,
                       ComplexObject<FunctionTree<1>> &inp,
                       int maxIter,
                       bool absPrec);

} // namespace mrcpp
