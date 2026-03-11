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

/**
 * @file HeatOperator.h
 * @brief Declaration of a separable convolution operator that realizes the heat
 *        semigroup \( e^{t\Delta} \) in \(D\) dimensions.
 *
 * @details
 * In \f$\mathbb{R}^D\f$, the heat propagator at time \(t>0\) is a Gaussian
 * convolution
 * \f[
 *   (e^{t\Delta} f)(\mathbf x)
 *   =
 *   \int_{\mathbb{R}^D}
 *     K_t(\mathbf x-\mathbf y)\, f(\mathbf y)\, d\mathbf y,
 *   \qquad
 *   K_t(\mathbf r)
 *   =
 *   \frac{1}{(4\pi t)^{D/2}}
 *   \exp\!\left(-\frac{\|\mathbf r\|^2}{4t}\right).
 * \f]
 * This class builds a rank-1 separable @ref ConvolutionOperator using a single 1D
 * Gaussian kernel in each coordinate and assembles the \(D\)-dimensional operator
 * as their tensor product. The amplitude/exponent are chosen so the overall kernel
 * matches \(K_t\).
 *
 * Construction delegates to the base class to:
 *  - project the 1D kernel to a function tree on a 1D MRA,
 *  - lift it to operator trees via cross-correlation,
 *  - transform/caches the result in the multiwavelet domain.
 *
 * The overload with explicit @p root and @p reach is useful for periodic boundary
 * conditions (PBC) or when the operator must be confined to a specific scale window.
 *
 * @see ConvolutionOperator, HeatKernel
 */

/**
 * @class HeatOperator
 * @ingroup operators
 * @brief D-dimensional heat semigroup as a separable Gaussian convolution.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 *
 * @note The kernel is normalized so that \(\int_{\mathbb{R}^D} K_t = 1\) and
 *       the map \(f \mapsto e^{t\Delta}f\) is positivity-preserving and
 *       \(L^1\)-contractive in the continuous setting.
 */
template <int D> class HeatOperator final : public ConvolutionOperator<D> {
public:
    /**
     * @brief Construct the heat operator \(e^{t\Delta}\) on the default scale window.
     *
     * @param mra  D-dimensional @ref MultiResolutionAnalysis defining domain and basis.
     * @param t    Diffusion time; must be strictly positive. Smaller @p t yields a
     *             narrower Gaussian and requires finer resolution.
     * @param prec Target build precision used while projecting the kernel and assembling
     *             the operator.
     *
     * @pre @p t > 0.
     * @see ConvolutionOperator
     */
    HeatOperator(const MultiResolutionAnalysis<D> &mra, double t, double prec);

    /**
     * @brief Construct the heat operator with an explicit operator scale window.
     *
     * @param mra   D-dimensional @ref MultiResolutionAnalysis.
     * @param t     Diffusion time; must be strictly positive.
     * @param prec  Target build precision.
     * @param root  Operator root (coarsest) scale.
     * @param reach Operator bandwidth (half-width in levels) at @p root; useful for
     *              periodic boundary conditions or domain tiling. Defaults to 1.
     *
     * @pre @p t > 0.
     * @see MWOperator, ConvolutionOperator
     */
    HeatOperator(const MultiResolutionAnalysis<D> &mra, double t, double prec, int root, int reach = 1);

    HeatOperator(const HeatOperator &oper) = delete;            ///< Non-copyable
    HeatOperator &operator=(const HeatOperator &oper) = delete; ///< Non-assignable
};

} // namespace mrcpp