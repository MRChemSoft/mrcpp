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

#include "ConvolutionOperator.h"

namespace mrcpp {

/** @class DerivativeConvolution
 * @ingroup operators
 *
 * @brief Separable convolution operator that approximates a spatial derivative
 *        using a differentiated Gaussian kernel.
 *
 * @tparam D Spatial dimension of the target operator (1, 2, or 3).
 *
 * @details
 * This class implements a *proof-of-concept* derivative as a convolution with the
 * derivative of a narrow Gaussian. In distributional terms one would like to have
 * \f$ \partial_x \delta \f$; numerically, we approximate it by
 * a derivative-of-Gaussian (DoG) kernel that is narrow enough to capture the
 * local derivative while remaining representable on the multiwavelet grid.
 *
 * Formally, for one Cartesian component:
 * \f[
 *   (D^x f)(\mathbf r)
 *   \;\approx\;
 *   \int_{\mathbb R^D}
 *     \frac{\partial}{\partial x}
 *     \left[\alpha\, e^{-\beta \lvert \mathbf r - \mathbf r' \rvert^2}\right]
 *     f(\mathbf r') \, d\mathbf r'
 *   \;=\; (k'_x * f)(\mathbf r).
 * \f]
 * In MRCPP this is realized as a @ref ConvolutionOperator with a 1D DoG kernel;
 * the D-dimensional operator is assembled as a separable tensor product across
 * coordinates and lifted to the multiwavelet basis.
 *
 * ### Responsibilities and division of labor
 * - **This class**: chooses a derivative-like kernel and exposes convenient
 *   constructors. It does not change application logic.
 * - **@ref ConvolutionOperator**: projects the 1D kernel to a function tree,
 *   lifts to operator trees via cross-correlation, transforms/caches in the MW basis,
 *   and manages rank/separability.
 *
 * ### Precision handling
 * The constructors accept a *build precision* that governs the narrowness of the
 * DoG kernel and the tolerances used during kernel projection and operator assembly:
 * tighter precision ⇒ narrower kernel ⇒ closer to an ideal derivative, but with
 * higher resolution demands and potentially larger operator bandwidth.
 *
 * @note This operator is primarily for validation/experiments. For production
 *       use consider:
 *       - @ref ABGVOperator for cuspy/discontinuous functions.
 *       - @ref BSOperator for sufficiently smooth functions.
 *
 * @see ConvolutionOperator, ABGVOperator, BSOperator
 */
template <int D> class DerivativeConvolution final : public ConvolutionOperator<D> {
public:
    /**
     * @brief Build a derivative-convolution operator on the default root/reach.
     *
     * @param mra  D-dimensional @ref MultiResolutionAnalysis that defines the domain
     *             and the scaling basis used by the operator.
     * @param prec Target build precision controlling kernel narrowness and
     *             assembly tolerances (tighter ⇒ narrower DoG).
     *
     * @details
     * Internally constructs a single-term derivative kernel (derivative of a Gaussian)
     * with parameters derived from @p prec, then delegates to
     * @ref ConvolutionOperator::initialize to:
     *  - project the 1D kernel to a function tree,
     *  - lift it to separable operator blocks via cross-correlation,
     *  - transform and cache in the multiwavelet basis.
     *
     * @warning Very small @p prec values produce *very* narrow kernels that may
     *          require deeper trees and higher-order bases to avoid under-resolution.
     */
    DerivativeConvolution(const MultiResolutionAnalysis<D> &mra, double prec);

    /**
     * @brief Build a derivative-convolution operator with an explicit scale window.
     *
     * @param mra   D-dimensional @ref MultiResolutionAnalysis.
     * @param prec  Target build precision (as above).
     * @param root  Operator root level (coarsest active scale).
     * @param reach Operator reach in levels; negative values trigger auto-detection
     *              from the domain extents.
     *
     * @details
     * Use this overload to constrain the operator to a chosen scale window; useful for
     * benchmarking, domain-decomposition experiments, or when composing multiple
     * operators with controlled supports. Kernel choice and assembly otherwise mirror
     * the simpler constructor.
     */
    DerivativeConvolution(const MultiResolutionAnalysis<D> &mra, double prec, int root, int reach);

    DerivativeConvolution(const DerivativeConvolution &oper) = delete;            ///< Non-copyable
    DerivativeConvolution &operator=(const DerivativeConvolution &oper) = delete; ///< Non-assignable
};

} // namespace mrcpp