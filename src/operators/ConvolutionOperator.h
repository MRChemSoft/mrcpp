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

#include "MWOperator.h"

namespace mrcpp {

/** @class ConvolutionOperator
 *  @ingroup operators
 *
 * @brief D-dimensional separable convolution operator built from a 1D Gaussian expansion.
 *
 * @tparam D Spatial dimension of the target operator.
 *
 * @details
 * This operator represents a separable convolution constructed from a sum of
 * one–dimensional Gaussian factors:
 *
 * \f[
 *    T \;=\; \sum_{m=1}^{M}
 *    \operatorname{sign}(\alpha_m)
 *    \bigotimes_{d=1}^{D}
 *    T_d\!\left(\beta_m,\;\sqrt[D]{|\alpha_m|}\right),
 * \f]
 *
 * where each \f$ T_d(\beta,\alpha) \f$ is the 1D convolution with kernel
 * \f$ k(x_d) = \alpha \exp(-\beta x_d^2) \f$ along coordinate \f$ x_d \f$.
 * The separable rank of the constructed operator equals the number of terms
 * \f$ M \f$ in the 1D Gaussian expansion:
 *
 * \f[
 *    \sum_{m=1}^{M} \alpha_m \exp(-\beta_m |x|^2).
 * \f]
 *
 * ### Construction strategy (high level)
 * 1. For each Gaussian term \f$ \alpha_m e^{-\beta_m x^2} \f$ in the 1D expansion,
 *    we rescale its coefficient to \f$ \sqrt[D]{|\alpha_m|} \f$ and keep the sign,
 *    so that the D-fold separable composition recovers the desired amplitude.
 * 2. Each 1D term is projected to a 1D multiresolution function tree using the
 *    same scaling family as the D-D operator (interpolating or Legendre).
 * 3. Cross-correlation machinery lifts each 1D factor into a 2D operator block
 *    (per axis-pair) and the @ref MWOperator backbone assembles the full D-D,
 *    separable operator.
 *
 * ### Precision control
 * The *build precision* (see @ref setBuildPrec and @ref getBuildPrec) governs:
 *  - the tolerance used when projecting 1D kernel terms to their function trees, and
 *  - the tolerance for assembling/thresholding operator trees.
 * Implementations typically use a tighter internal precision for the kernel
 * projection than for the operator assembly to keep the total error within the
 * requested target.
 *
 * @note All constructors are *non-owning* with respect to the input expansion; the
 * implementation copies kernel terms as needed for projection, and keeps only the
 * operator trees internally.
 *
 * @see ConvolutionOperator::initialize
 * @see ConvolutionOperator::getKernelMRA
 * @see MWOperator
 */
template <int D> class ConvolutionOperator : public MWOperator<D> {
public:
    /**
     * @brief Build a separable convolution operator on the default operator root/extent.
     *
     * @param mra     D-dimensional @ref MultiResolutionAnalysis defining the domain and basis.
     * @param kernel  1D Gaussian expansion providing the separable factors (rank = kernel.size()).
     * @param prec    Target build precision used to steer kernel projection and operator assembly.
     *
     * @details Uses the operator's default root scale (@c mra.getRootScale()) and a
     * reach chosen by the implementation. For more control over root/reach, use the
     * other constructor.
     */
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, GaussExp<1> &kernel, double prec);

    /**
     * @brief Build a separable convolution operator with explicit root scale and reach.
     *
     * @param mra     D-dimensional @ref MultiResolutionAnalysis.
     * @param kernel  1D Gaussian expansion (rank = kernel.size()).
     * @param prec    Target build precision.
     * @param root    Operator root scale (level) to anchor the construction.
     * @param reach   Operator reach (half-width in levels). Negative value means
     *                *auto*—deduced from the world box extents.
     *
     * @details This variant allows advanced users to control the scale window spanned
     * by the operator representation, which may be useful when coupling to other
     * operators or enforcing boundary extents.
     */
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, GaussExp<1> &kernel, double prec, int root, int reach);

    ConvolutionOperator(const ConvolutionOperator &oper) = delete;
    ConvolutionOperator &operator=(const ConvolutionOperator &oper) = delete;
    virtual ~ConvolutionOperator() = default;

    /// @brief Retrieve the user-requested build precision associated with this operator.
    double getBuildPrec() const { return this->build_prec; }

protected:
    /**
     * @brief Protected convenience constructor for subclasses that defer initialization.
     *
     * @param mra D-dimensional @ref MultiResolutionAnalysis.
     *
     * @details Initializes the @ref MWOperator base with default root and reach.
     * Subclasses must call @ref initialize to populate the separable expansion.
     */
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra)
            : MWOperator<D>(mra, mra.getRootScale(), -10) {}

    /**
     * @brief Protected convenience constructor with explicit root and reach.
     *
     * @param mra   D-dimensional @ref MultiResolutionAnalysis.
     * @param root  Root level.
     * @param reach Reach (half-width in levels). Negative = auto.
     */
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, int root, int reach)
            : MWOperator<D>(mra, root, reach) {}

    /**
     * @brief Core build routine that projects 1D kernel terms and assembles operator trees.
     *
     * @param kernel 1D Gaussian expansion (input rank M).
     * @param k_prec Precision used when projecting each 1D Gaussian term into a
     *               1D function tree (typically tighter than @p o_prec).
     * @param o_prec Precision used when expanding to operator trees and performing
     *               wavelet transforms/thresholding.
     *
     * @details For each term in @p kernel:
     *  - Coefficient is rescaled to \f$ \sqrt[D]{|\alpha|} \f$ with the original sign,
     *    ensuring the D-fold separable product reproduces the intended amplitude.
     *  - The analytic 1D Gaussian is projected to a 1D @ref FunctionTree with tolerance
     *    @p k_prec.
     *  - A @ref CrossCorrelationCalculator lifts the 1D representation to a 2D operator
     *    block; bottom-up wavelet transforms and caching finalize each block.
     * The set of blocks is stored in the @ref MWOperator base and exposed as a
     * separable expansion of rank @c kernel.size().
     */
    void initialize(GaussExp<1> &kernel, double k_prec, double o_prec);

    /// @brief Store the user-requested build precision (used for reporting/inspection).
    void setBuildPrec(double prec) { this->build_prec = prec; }

    /**
     * @brief Build a 1D @ref MultiResolutionAnalysis to discretize kernel factors.
     *
     * @return A 1D MRA whose scaling family matches the D-D operator MRA (Interpolating or Legendre),
     *         with an order chosen as \f$ 2s+1 \f$ where \f$ s \f$ is the operator scaling order.
     *
     * @details The 1D box uses the operator root. Its reach is the operator reach + 1
     * (or derived from the world box if reach is negative) to ensure kernel support
     * covers the correlations used during lifting.
     */
    MultiResolutionAnalysis<1> getKernelMRA() const;

    /// Target precision requested at construction time; used to steer sub-steps in the build.
    double build_prec{-1.0};
};

} // namespace mrcpp
