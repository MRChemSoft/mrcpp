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
 * @file IdentityConvolution.cpp
 * @brief Implementation of a separable convolution operator that approximates the identity
 *        via a narrow Gaussian kernel (delta-approximation).
 *
 * @details
 * This translation unit defines the templated constructors of
 * @ref mrcpp::IdentityConvolution, a convenience @ref ConvolutionOperator that uses a
 * single-term Gaussian kernel to approximate the Dirac delta distribution:
 * \f[
 *   \delta(x) \;\approx\; \alpha\,e^{-\beta x^2}.
 * \f]
 * The associated D-dimensional operator is assembled separably (tensor-product form)
 * following MRCPP’s multiwavelet machinery. The build precision controls the kernel
 * narrowness and the tolerances used during projection/assembly.
 */

#include "IdentityConvolution.h"
#include "IdentityKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @brief Construct an identity-like convolution operator on the default scale window.
 *
 * @tparam D   Spatial dimension (1, 2, or 3).
 * @param mra  D-dimensional @ref MultiResolutionAnalysis that defines the domain and basis.
 * @param prec Target build precision controlling the closeness to the delta function.
 *
 * @details
 * Internally the constructor:
 *  - Stores @p prec as the build precision.
 *  - Uses split tolerances:
 *      - @c k_prec = prec/10.0 for accurate projection of the narrow Gaussian kernel.
 *      - @c o_prec = prec for operator assembly.
 *  - Builds a single-term @ref IdentityKernel<D> at @c k_prec and calls
 *    @ref ConvolutionOperator::initialize to lift it into separable operator blocks.
 *  - Finalizes with @ref MWOperator::initOperExp for bookkeeping/caching.
 *
 * A tighter @p prec yields a narrower Gaussian (better delta approximation) but
 * increases the required resolution and operator bandwidth in practice.
 */
template <int D>
IdentityConvolution<D>::IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double prec)
        : ConvolutionOperator<D>(mra) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 10.0;

    IdentityKernel<D> kernel(k_prec);
    this->initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Construct an identity-like convolution operator with explicit root and reach (PBC-ready).
 *
 * @tparam D    Spatial dimension (1, 2, or 3).
 * @param mra   D-dimensional @ref MultiResolutionAnalysis.
 * @param prec  Target build precision controlling the closeness to the delta function.
 * @param root  Operator root level (coarsest scale at which the operator is defined).
 * @param reach Operator half-bandwidth at @p root (relevant for periodic boundary conditions).
 *
 * @details
 * This overload confines the operator to a specific scale window, which is useful for
 * periodic boundary conditions or when coupling multiple operators with controlled support.
 * Compared to the default constructor, the kernel projection tolerance is chosen even
 * tighter (@c k_prec = prec/100.0) to ensure faithful representation on restricted scale
 * ranges; operator assembly uses @c o_prec = prec.
 *
 * Steps:
 *  1. Record @p prec via @ref ConvolutionOperator::setBuildPrec.
 *  2. Create a single-term @ref IdentityKernel<D> at @c k_prec.
 *  3. Initialize separable operator blocks (@ref ConvolutionOperator::initialize)
 *     within the user-specified scale window (@p root, @p reach).
 *  4. Call @ref MWOperator::initOperExp.
 */
template <int D>
IdentityConvolution<D>::IdentityConvolution(const MultiResolutionAnalysis<D> &mra, double prec, int root, int reach)
        : ConvolutionOperator<D>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    double o_prec = prec;
    double k_prec = prec / 100.0;

    IdentityKernel<D> kernel(k_prec);
    this->initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

/* Explicit template instantiations */
template class IdentityConvolution<1>;
template class IdentityConvolution<2>;
template class IdentityConvolution<3>;

} // namespace mrcpp