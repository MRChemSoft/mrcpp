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
 * @file PoissonOperator.h
 * @brief Separable multiwavelet convolution operator for the 3D Poisson kernel.
 *
 * The operator realizes a fast approximation of the Green's function
 * \f$ P(\mathbf r-\mathbf r') = 1/\lvert \mathbf r-\mathbf r'\rvert \f$
 * by expanding it into a finite sum of Gaussians,
 * \f[
 *   \frac{1}{\lvert \mathbf r-\mathbf r'\rvert}
 *   \;\approx\;
 *   \sum_{m=1}^{M} \alpha_m \exp\!\big(-\beta_m \lvert \mathbf r-\mathbf r'\rvert^2\big),
 * \f]
 * which enables a tensor–separable application along Cartesian axes in the
 * multiwavelet framework. See @ref ConvolutionOperator for assembly details and
 * @ref PoissonOperator (implementation) for construction mechanics.
 */

#pragma once

#include "ConvolutionOperator.h"

namespace mrcpp {

/**
 * @class PoissonOperator
 * @brief Convolution with the Poisson Green's function kernel in 3D.
 *
 * @details
 * The Poisson kernel is approximated by a Gaussian expansion, allowing the operator
 * to be applied as a separated product over Cartesian directions:
 * \f[
 *   P(\mathbf r-\mathbf r')
 *   = \frac{1}{\lvert \mathbf r-\mathbf r'\rvert}
 *   \;\approx\; \sum_{m=1}^{M} \alpha_m \exp\!\big(-\beta_m \lvert \mathbf r-\mathbf r'\rvert^2\big).
 * \f]
 * Each 1D Gaussian term is projected to a function tree and lifted to operator blocks
 * via cross-correlation; the full 3D operator is then cached for efficient application.
 *
 * The expansion accuracy and kernel width are controlled by the requested build precision.
 * An overload with explicit @p root/@p reach can confine the operator to a chosen scale window
 * (useful for periodic-style setups or domain-decomposition experiments).
 *
 * @see ConvolutionOperator, PoissonKernel
 */
class PoissonOperator final : public ConvolutionOperator<3> {
public:
    /**
     * @brief Construct a Poisson operator on the default root/reach of the provided MRA.
     *
     * @param mra   3D @ref MultiResolutionAnalysis defining the domain and scaling basis.
     * @param prec  Target build precision controlling the Gaussian expansion (smaller ⇒ tighter/longer rank).
     */
    PoissonOperator(const MultiResolutionAnalysis<3> &mra, double prec);

    /**
     * @brief Construct a Poisson operator with an explicit scale window.
     *
     * @param mra    3D @ref MultiResolutionAnalysis.
     * @param prec   Target build precision.
     * @param root   Operator root level (coarsest scale where the operator resides).
     * @param reach  Operator reach (half-width, in levels) around @p root; affects bandwidth/PBC-like extent.
     */
    PoissonOperator(const MultiResolutionAnalysis<3> &mra, double prec, int root, int reach = 1);

    PoissonOperator(const PoissonOperator &oper) = delete;            ///< Non-copyable
    PoissonOperator &operator=(const PoissonOperator &oper) = delete; ///< Non-assignable
};

} // namespace mrcpp