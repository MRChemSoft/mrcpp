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
 * @file HeatOperator.cpp
 * @brief Implementation of a separable convolution operator that realizes the
 *        \(D\)-dimensional heat semigroup via a single-term Gaussian kernel.
 *
 * @details
 * The continuous heat propagator at time \(t>0\) is the convolution with
 * \f[
 *   K_t(\mathbf x)
 *   =
 *   \frac{1}{(4\pi t)^{D/2}}
 *   \exp\!\left(-\frac{\lVert\mathbf x\rVert^2}{4t}\right).
 * \f]
 * In MRCPP, separable operators are assembled from 1D Gaussian building blocks.
 * This implementation:
 *  - constructs a @ref HeatKernel whose single 1D Gaussian has exponent
 *    \f$\beta = 1/(4t)\f$ and coefficient \f$\alpha = (\beta/\pi)^{D/2}\f$,
 *  - projects that kernel to a 1D function tree,
 *  - lifts it to an operator tree by cross-correlation,
 *  - transforms/caches the operator in the multiwavelet domain,
 *  - and exposes it as a @ref ConvolutionOperator acting in \(D\) dimensions.
 *
 * Two constructors are provided: a default one (using the operator's default
 * root/reach) and one tailored for periodic boundary conditions (PBC) with an
 * explicit scale window @p root/@p reach.
 */

#include "HeatOperator.h"
#include "HeatKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @brief Build a heat propagator \(e^{t\Delta}\) as a separable convolution.
 *
 * @tparam D Spatial dimension of the operator (1, 2, or 3).
 *
 * @param[in] mra  D-dimensional @ref MultiResolutionAnalysis that defines both
 *                 the computational domain and the scaling basis.
 * @param[in] t    Diffusion time; determines the kernel width
 *                 (\f$\beta = 1/(4t)\f$). Must be strictly positive.
 * @param[in] prec Target build precision for assembling the operator.
 *
 * @details
 * Steps performed:
 *  1. The requested build precision is recorded via @c setBuildPrec(prec).
 *  2. Two tolerances are chosen:
 *     - @c k_prec = @p prec / 10 for the 1D kernel projection (tighter),
 *     - @c o_prec = @p prec for the operator assembly.
 *  3. A @ref HeatKernel<D> is instantiated with exponent \f$1/(4t)\f$ and
 *     amplitude chosen to match \f$(4\pi t)^{-D/2}\f$ upon separable assembly.
 *  4. @ref ConvolutionOperator::initialize is called to:
 *     - project the kernel to a 1D function tree,
 *     - lift it to an operator tree via cross-correlation,
 *     - transform and cache the operator in the MW domain.
 *  5. @ref initOperExp is called to finalize the separable expansion (rank 1).
 *
 * @note Smaller @p t \(\Rightarrow\) narrower Gaussian \(\Rightarrow\) more demanding
 *       resolution (consider tightening @p prec and/or extending operator reach).
 * @warning Passing non-positive @p t yields a meaningless kernel; callers must
 *          ensure @p t > 0.
 */
template <int D>
HeatOperator<D>::HeatOperator(const MultiResolutionAnalysis<D> &mra, double t, double prec)
        : ConvolutionOperator<D>(mra) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;        // Operator-assembly tolerance
    double k_prec = prec / 10.0; // Kernel-projection tolerance (tighter)

    HeatKernel<D> kernel(t);
    this->initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size()); // single-term expansion (rank = 1)

    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Build a heat propagator with an explicit operator scale window (PBC use).
 *
 * @tparam D Spatial dimension of the operator (1, 2, or 3).
 *
 * @param[in] mra   D-dimensional @ref MultiResolutionAnalysis.
 * @param[in] t     Diffusion time (\f$t>0\f$).
 * @param[in] prec  Target build precision.
 * @param[in] root  Root (coarsest) scale the operator is attached to.
 * @param[in] reach Bandwidth at the root scale (useful for PBC/domain tiling).
 *
 * @details
 * This overload mirrors the default constructor but confines the operator to a
 * specific scale window, which is particularly useful for periodic boundary
 * conditions and domain-decomposition setups.
 *
 * Implementation differences vs. the default constructor:
 *  - The base @ref ConvolutionOperator is constructed with (@p root, @p reach).
 *  - @c k_prec is chosen even tighter ( @p prec / 100.0 ) to robustly capture
 *    the narrow Gaussian under potentially coarser scale constraints.
 *
 * @note The @p reach parameter controls the operator bandwidth measured in
 *       levels at @p root; see @ref MWOperator for details on scale windows
 *       and bandwidth semantics.
 */
template <int D>
HeatOperator<D>::HeatOperator(const MultiResolutionAnalysis<D> &mra, double t, double prec, int root, int reach)
        : ConvolutionOperator<D>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;          // Operator-assembly tolerance
    double k_prec = prec / 100.0;  // Very tight kernel-projection tolerance

    HeatKernel<D> kernel(t);
    this->initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size()); // single-term expansion (rank = 1)

    Printer::setPrintLevel(oldlevel);
}

/* Explicit template instantiations */
template class HeatOperator<1>;
template class HeatOperator<2>;
template class HeatOperator<3>;

} // namespace mrcpp