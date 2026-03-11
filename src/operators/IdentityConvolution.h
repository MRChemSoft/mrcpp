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
 * @file IdentityConvolution.h
 * @brief Separable convolution operator that approximates the identity (Dirac delta)
 *        using a narrow Gaussian kernel.
 *
 * @details
 * This header declares @ref mrcpp::IdentityConvolution, a thin convenience wrapper
 * around @ref mrcpp::ConvolutionOperator that realizes an identity-like operator via
 * a single Gaussian kernel, separably assembled in D dimensions.
 *
 * The kernel approximation is
 * \f[
 *   \delta(\mathbf r - \mathbf r')
 *   \;\approx\;
 *   \alpha \exp\!\bigl(-\beta\lVert \mathbf r - \mathbf r' \rVert^2\bigr),
 * \f]
 * which, in the MRCPP framework, is projected to a 1D function tree and lifted to
 * D-dimensional operator blocks by cross-correlation. The resulting operator is
 * bandwidth-limited and numerically stable for use in multiresolution workflows.
 *
 * The constructor takes a *build precision* that governs the kernel’s narrowness
 * and the tolerances used during projection and assembly. Tighter precision yields
 * a Gaussian closer to a true delta (hence a better identity approximation), at the
 * cost of higher resolution demands.
 */

#pragma once

#include "ConvolutionOperator.h"

namespace mrcpp {

/**
 * @class IdentityConvolution
 * @ingroup operators
 * @brief Convolution with an identity (delta-like) kernel.
 *
 * @tparam D Spatial dimension of the target operator (1, 2, or 3).
 *
 * @details
 * The operator is represented as a separable sum (rank-1 in the default realization)
 * of 1D Gaussian convolutions identical along each Cartesian direction. It is mainly
 * intended for diagnostics and algorithmic baselines; for strict identity action,
 * prefer direct coefficient transfers when applicable.
 *
 * The underlying kernel is the Gaussian surrogate of the Dirac delta,
 * \f$ I(\mathbf r-\mathbf r') \approx \alpha e^{-\beta \lVert \mathbf r-\mathbf r' \rVert^2} \f$,
 * with parameters chosen from the requested build precision.
 *
 * @see ConvolutionOperator
 */
template <int D> class IdentityConvolution final : public ConvolutionOperator<D> {
public:
    /**
     * @brief Build an identity-like convolution operator on the default root/reach.
     *
     * @param mra  D-dimensional @ref MultiResolutionAnalysis defining domain and basis.
     * @param prec Target build precision controlling the closeness to the delta function
     *             (narrowness of the Gaussian) and assembly tolerances.
     *
     * @details
     * Internally constructs a single-term Gaussian kernel and invokes
     * @ref ConvolutionOperator::initialize to assemble the separable operator blocks.
     */
    IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double prec);

    /**
     * @brief Build an identity-like convolution operator with explicit scale window.
     *
     * @param mra   D-dimensional @ref MultiResolutionAnalysis.
     * @param prec  Target build precision (as above).
     * @param root  Operator root level (coarsest scale at which the operator resides).
     * @param reach Operator half-bandwidth at @p root (useful for periodic domains).
     *
     * @details
     * Use this overload to confine the operator to a specific scale window—particularly
     * helpful under periodic boundary conditions or when coordinating multiple operators.
     */
    IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double prec, int root, int reach = 1);

    IdentityConvolution(const IdentityConvolution &oper) = delete;            ///< Non-copyable
    IdentityConvolution &operator=(const IdentityConvolution &oper) = delete; ///< Non-assignable
};

} // namespace mrcpp