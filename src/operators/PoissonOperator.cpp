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
 * @file PoissonOperator.cpp
 * @brief Definition of a separable convolution operator that approximates the 3D Poisson kernel.
 *
 * The operator is assembled from a one–dimensional Gaussian expansion of \f$1/r\f$
 * (see @ref PoissonKernel). Each 1D term is projected to a function tree and lifted
 * to a 2D operator block by cross-correlation; the full 3D operator is built as a
 * separable product and cached for efficient application.
 */

#include "PoissonOperator.h"
#include "PoissonKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @class PoissonOperator
 * @brief Convolution operator approximating the 3D Poisson kernel \f$1/\lvert \mathbf r\rvert\f$.
 *
 * @details
 * The kernel is approximated on a bounded radial interval by a finite Gaussian expansion
 * \f[
 *   \frac{1}{r} \approx \sum_{m=1}^{M} \beta_m\, e^{-\alpha_m r^2},
 * \f]
 * which enables a separated representation amenable to fast multiwavelet application
 * along Cartesian axes. Construction proceeds by:
 * 1) choosing a target build precision to set the effective kernel width,
 * 2) computing a validity interval \f$[r_{\min}, r_{\max}]\f$ from the MRA,
 * 3) generating the Gaussian terms via @ref PoissonKernel, and
 * 4) projecting and lifting each term into operator blocks before caching.
 */

/**
 * @brief Build a Poisson operator on the default root/reach of the provided MRA.
 *
 * @param mra   Three–dimensional @ref MultiResolutionAnalysis defining domain and basis.
 * @param prec  Target build precision (heuristic closeness to \f$1/r\f$); smaller ⇒ tighter kernel.
 *
 * @details
 * - Uses @c k_prec = prec/10 for kernel projection and @c o_prec = prec for operator assembly.
 * - The radial interval is inferred from @p mra:
 *   - \f$ r_{\min} = \text{MRA.calcMinDistance}(k\_prec) \f$ (resolution-limited),
 *   - \f$ r_{\max} = \text{MRA.calcMaxDistance}() \f$ (domain-limited).
 * - Constructs a @ref PoissonKernel on \f$[r_{\min}, r_{\max}]\f$, initializes internal
 *   operator trees, and prepares caches for application.
 */
PoissonOperator::PoissonOperator(const MultiResolutionAnalysis<3> &mra, double prec)
        : ConvolutionOperator<3>(mra) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;         // operator-assembly tolerance
    double k_prec = prec / 10.0;  // kernel-projection tolerance
    double r_min = this->MRA.calcMinDistance(k_prec);
    double r_max = this->MRA.calcMaxDistance();

    PoissonKernel kernel(k_prec, r_min, r_max);
    initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Build a Poisson operator with explicit scale window (root/reach), e.g. for PBC-style setups.
 *
 * @param mra    Three–dimensional @ref MultiResolutionAnalysis.
 * @param prec   Target build precision (heuristic closeness to \f$1/r\f$).
 * @param root   Operator root level (coarsest scale where the operator lives).
 * @param reach  Operator reach in levels (half-width around @p root); negative ⇒ auto-detect.
 *
 * @details
 * - Uses a tighter kernel projection tolerance @c k_prec = prec/100 and @c o_prec = prec
 *   for assembling the operator blocks.
 * - The base radial extent is obtained from the MRA; then \f$ r_{\max} \f$ is rescaled
 *   to reflect the selected operator scale window (periodic-world style adjustment):
 *   \f[
 *     r_{\max} \leftarrow r_{\max} \, 2^{-(\text{oper\_root} - \text{MRA.getRootScale()})}
 *                         \, \bigl( 2\,\text{oper\_reach} + 1 \bigr).
 *   \f]
 * - Constructs and initializes the Gaussian expansion accordingly and prepares
 *   the operator components and caches.
 */
PoissonOperator::PoissonOperator(const MultiResolutionAnalysis<3> &mra, double prec, int root, int reach)
        : ConvolutionOperator<3>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;           // operator-assembly tolerance
    double k_prec = prec / 100.0;   // very tight kernel-projection tolerance
    double r_min = this->MRA.calcMinDistance(k_prec);
    double r_max = this->MRA.calcMaxDistance();

    // Adjust r_max to the chosen operator scale window (e.g., periodic-world bandwidth)
    auto rel_root = this->oper_root - this->MRA.getRootScale();
    r_max *= std::pow(2.0, -rel_root);
    r_max *= (2.0 * this->oper_reach) + 1.0;

    PoissonKernel kernel(k_prec, r_min, r_max);
    initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

} // namespace mrcpp