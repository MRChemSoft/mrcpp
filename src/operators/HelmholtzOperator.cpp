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
 * @file HelmholtzOperator.cpp
 * @brief Definition of a separable 3D convolution operator approximating the Helmholtz/Yukawa kernel.
 *
 * @details
 * This file implements @ref mrcpp::HelmholtzOperator, a convenience convolution operator
 * in three spatial dimensions that applies a Gaussian expansion of the radial kernel
 * \f$ e^{-\mu r}/r \f$. The expansion is built by @ref mrcpp::HelmholtzKernel and
 * lifted into a separable multiwavelet operator, which can then be applied along the
 * Cartesian directions.
 */

#include "HelmholtzOperator.h"
#include "HelmholtzKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @brief Construct a 3D Helmholtz (Yukawa) convolution operator.
 *
 * @param mra   The 3D @ref MultiResolutionAnalysis that defines the domain and basis.
 * @param mu    Screening parameter \f$\mu>0\f$ of the Yukawa kernel.
 * @param prec  Build precision; controls kernel- and operator-assembly tolerances and,
 *              indirectly, the separation rank of the Gaussian expansion.
 *
 * @details
 * - Chooses separate tolerances for the kernel projection (@c k_prec = prec/10) and
 *   the operator assembly (@c o_prec = prec).
 * - Estimates the admissible radial interval \f$[r_{\min}, r_{\max}]\f$ from @p mra via
 *   @ref MultiResolutionAnalysis::calcMinDistance and @ref MultiResolutionAnalysis::calcMaxDistance.
 * - Builds a @ref HelmholtzKernel on that interval with the requested accuracy, then
 *   calls @ref ConvolutionOperator::initialize to form the separable operator blocks
 *   and caches them via @ref MWOperator::initOperExp.
 *
 * @note The printer level is temporarily reduced during build to keep output concise.
 */
HelmholtzOperator::HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double mu, double prec)
        : ConvolutionOperator<3>(mra) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 10.0;
    double r_min = this->MRA.calcMinDistance(k_prec);
    double r_max = this->MRA.calcMaxDistance();

    HelmholtzKernel kernel(mu, k_prec, r_min, r_max);
    initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Construct a 3D Helmholtz convolution operator with explicit root and reach.
 *
 * @param mra    The 3D @ref MultiResolutionAnalysis.
 * @param mu     Screening parameter \f$\mu>0\f$.
 * @param prec   Build precision (as above).
 * @param root   Operator root scale (coarsest level for the operator support).
 * @param reach  Operator reach (half-width in levels). For periodic domains this
 *               sets the operator bandwidth at @p root.
 *
 * @details
 * - Uses a tighter kernel-projection tolerance (@c k_prec = prec/100) while keeping
 *   the operator-assembly tolerance at @c o_prec = prec.
 * - Estimates \f$[r_{\min}, r_{\max}]\f$ as in the other constructor, then adjusts
 *   @c r_max to reflect periodic worlds by scaling with the relative root shift and
 *   the chosen @p reach:
 *   \f[
 *     r_{\max} \leftarrow r_{\max}\, 2^{-(\text{oper\_root} - \text{MRA.root})}
 *               \times \big( 2\,\text{reach} + 1 \big).
 *   \f]
 * - Builds the @ref HelmholtzKernel and initializes the separable operator.
 *
 * @note This overload is intended for periodic boundary conditions or scenarios
 *       where the operator must be confined to a specific scale window.
 */
HelmholtzOperator::HelmholtzOperator(const MultiResolutionAnalysis<3> &mra, double mu, double prec, int root, int reach)
        : ConvolutionOperator<3>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 100.0;
    double r_min = this->MRA.calcMinDistance(k_prec);
    double r_max = this->MRA.calcMaxDistance();

    // Adjust r_max for periodic world
    auto rel_root = this->oper_root - this->MRA.getRootScale();
    r_max *= std::pow(2.0, -rel_root);
    r_max *= (2.0 * this->oper_reach) + 1.0;

    HelmholtzKernel kernel(mu, k_prec, r_min, r_max);
    initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

} // namespace mrcpp