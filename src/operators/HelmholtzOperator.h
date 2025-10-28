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
 * @file HelmholtzOperator.h
 * @brief Declaration of a 3D separable convolution operator for the Helmholtz/Yukawa kernel.
 *
 * @details
 * This header declares @ref mrcpp::HelmholtzOperator, a specialized
 * @ref ConvolutionOperator that applies the screened Coulomb (Yukawa) Green's function
 * in three spatial dimensions via a Gaussian expansion. The radial kernel
 * \f$ e^{-\mu r}/r \f$ is approximated as a finite sum of 1D Gaussians, enabling
 * separated application across Cartesian coordinates in the MRCPP multiwavelet basis.
 */

#pragma once

#include "ConvolutionOperator.h"

namespace mrcpp {

/**
 * @class HelmholtzOperator
 * @ingroup operators
 *
 * @brief Separable 3D convolution approximating the Helmholtz (Yukawa) Green's function.
 *
 * @details
 * The continuous kernel
 * \f[
 *   H(\mathbf r - \mathbf r') = \frac{e^{-\mu \lvert \mathbf r - \mathbf r' \rvert}}
 *                                     {\lvert \mathbf r - \mathbf r' \rvert}
 * \f]
 * is approximated by a Gaussian sum
 * \f[
 *   H(\mathbf r - \mathbf r')
 *   \;\approx\;
 *   \sum_{m=1}^{M} \alpha_m \exp\!\big( -\beta_m \lvert \mathbf r - \mathbf r' \rvert^2 \big),
 * \f]
 * which admits a *separable* representation in Cartesian coordinates, allowing the
 * operator to be assembled as a tensor product of 1D convolution blocks within the
 * MRCPP framework. The expansion coefficients \f$ \alpha_m, \beta_m \f$ and the
 * separation rank \f$ M \f$ are chosen internally based on the requested build
 * precision.
 *
 * ### Usage notes
 * - This class is a convenience wrapper that constructs the Gaussian expansion and
 *   the corresponding multiwavelet operator trees; application is handled by the
 *   @ref ConvolutionOperator base.
 * - For periodic worlds and explicit scale control, use the constructor that accepts
 *   @p root and @p reach (see below).
 *
 * @see ConvolutionOperator, HelmholtzKernel
 */
class HelmholtzOperator final : public ConvolutionOperator<3> {
public:
    /**
     * @brief Build a Helmholtz (Yukawa) convolution operator on the default scale window.
     *
     * @param mra  3D @ref MultiResolutionAnalysis defining domain and basis.
     * @param m    Screening parameter \f$ \mu > 0 \f$ of the Yukawa kernel.
     * @param prec Target build precision controlling the Gaussian expansion accuracy
     *             and operator assembly tolerances.
     *
     * @details
     * Internally:
     * 1. Estimates admissible radial bounds from @p mra.
     * 2. Constructs a Gaussian expansion for \f$ e^{-\mu r}/r \f$ at the requested accuracy.
     * 3. Lifts the 1D kernels to separable operator blocks and caches them.
     */
    HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double m, double prec);

    /**
     * @brief Build a Helmholtz operator with explicit root scale and reach (useful for PBC).
     *
     * @param mra   3D @ref MultiResolutionAnalysis.
     * @param m     Screening parameter \f$ \mu \f$ (same as above).
     * @param prec  Target build precision.
     * @param root  Operator root level (coarsest scale at which the operator is defined).
     * @param reach Operator half-bandwidth at @p root (controls extent; relevant for periodic worlds).
     *
     * @details
     * This overload confines the operator to a specified scale window and adjusts the
     * radial extent accordingly—suitable for periodic boundary conditions and scenarios
     * requiring strict bandwidth control.
     */
    HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double m, double prec, int root, int reach = 1);

    HelmholtzOperator(const HelmholtzOperator &oper) = delete;            ///< Non-copyable
    HelmholtzOperator &operator=(const HelmholtzOperator &oper) = delete; ///< Non-assignable
};

} // namespace mrcpp