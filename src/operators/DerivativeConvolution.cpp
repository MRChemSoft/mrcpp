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

#include "DerivativeConvolution.h"
#include "DerivativeKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @class DerivativeConvolution
 * @brief Separable convolution operator that approximates a (first) derivative
 *        using a differentiated Gaussian kernel.
 *
 * @tparam D Spatial dimension of the target operator.
 *
 * @details
 * This implementation is a thin wrapper around @ref ConvolutionOperator that:
 *  - **Chooses the kernel**: a *single* derivative-of-Gaussian (DoG) term whose
 *    width is set by a requested build precision.
 *  - **Projects the 1D kernel** to a function tree at a tight tolerance.
 *  - **Lifts** the 1D kernel into a D-dimensional operator via cross-correlation
 *    (separable assembly) and prepares it for application in the multiwavelet basis.
 *
 * The resulting operator is bandwidth-limited and numerically stable, offering a
 * controlled approximation to a spatial derivative. It is mainly for validation
 * and experimentation; for production derivatives consider @ref ABGVOperator
 * (cusps/discontinuities) or @ref BSOperator (smooth functions).
 */

/**
 * @brief Construct a derivative-convolution operator on the default root/reach.
 *
 * @param mra  D-dimensional @ref MultiResolutionAnalysis defining basis/domain.
 * @param prec Target build precision that controls kernel width and assembly thresholds.
 *
 * @details
 * Steps performed here:
 *  1. **Silence verbose output** during operator build by temporarily lowering the
 *     global print level (restored upon exit).
 *  2. **Record build precision** via @c setBuildPrec(prec). This precision is later
 *     available from @ref ConvolutionOperator::getBuildPrec.
 *  3. **Split tolerances** into:
 *     - @c k_prec = prec/10.0 for *kernel projection* (tighter; DoG is narrow),
 *     - @c o_prec = prec       for *operator assembly* (adequate once kernel is accurate).
 *  4. **Create the kernel**: a single-term @ref DerivativeKernel<D> parametrized by
 *     @c k_prec, which internally chooses the Gaussian exponent consistent with the
 *     requested accuracy.
 *  5. **Assemble the operator** by calling @ref ConvolutionOperator::initialize,
 *     which projects the kernel to a 1D function tree, lifts it to operator trees
 *     via cross-correlation, transforms to the MW domain, and caches nodes.
 *
 * The operator rank equals the number of 1D kernel terms; for this kernel it is 1.
 */
template <int D>
DerivativeConvolution<D>::DerivativeConvolution(const MultiResolutionAnalysis<D> &mra, double prec)
        : ConvolutionOperator<D>(mra) {
    // Keep the build quiet; restore the previous level at the end of scope.
    int oldlevel = Printer::setPrintLevel(0);

    // Store build precision on the base class for later diagnostics/inspection.
    this->setBuildPrec(prec);

    // Operator-assembly tolerance: used while expanding/lifting the kernel to operator trees.
    double o_prec = prec;

    // Kernel-projection tolerance: tighter than operator assembly to resolve a narrow DoG.
    double k_prec = prec / 10.0;

    // A single differentiated Gaussian tuned by k_prec.
    DerivativeKernel<D> kernel(k_prec);

    // Build the separable operator blocks from the 1D kernel.
    this->initialize(kernel, k_prec, o_prec);

    // Restore previous print level.
    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Construct a derivative-convolution operator with explicit root and reach.
 *
 * @param mra    D-dimensional @ref MultiResolutionAnalysis.
 * @param prec   Target build precision that controls kernel width and thresholds.
 * @param root   Operator root level (coarsest active scale).
 * @param reach  Operator reach (half-width in levels). Negative => auto-detected.
 *
 * @details
 * This overload is identical in spirit to the simpler constructor, but allows
 * **explicit control of the active scale window**:
 *  - Use when benchmarking, debugging, or composing multiple operators whose
 *    supports must be constrained not to overlap.
 *
 * Implementation notes mirror the first ctor with one change:
 *  - The kernel projection is made *even tighter*: @c k_prec = prec/100.0,
 *    which helps when the operator is restricted to a narrower scale window
 *    (ensuring the DoG is still faithfully represented).
 */
template <int D>
DerivativeConvolution<D>::DerivativeConvolution(const MultiResolutionAnalysis<D> &mra, double prec, int root, int reach)
        : ConvolutionOperator<D>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);

    // Operator-assembly tolerance (same rationale as the other ctor).
    double o_prec = prec;

    // Very tight kernel-projection tolerance for explicit windowing.
    double k_prec = prec / 100.0;

    DerivativeKernel<D> kernel(k_prec);
    this->initialize(kernel, k_prec, o_prec);

    Printer::setPrintLevel(oldlevel);
}

/* Explicit template instantiations */
template class DerivativeConvolution<1>;
template class DerivativeConvolution<2>;
template class DerivativeConvolution<3>;

} // namespace mrcpp