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

#include "functions/GaussExp.h"
#include "functions/GaussFunc.h"
#include "functions/GaussPoly.h"

namespace mrcpp {

/**
 * @class DerivativeKernel
 * @ingroup operators
 *
 * @brief One–dimensional *derivative-of-Gaussian* (DoG) kernel packaged as a
 *        Gaussian expansion of rank 1, suitable for building separable
 *        convolution-based derivative operators in D dimensions.
 *
 * @tparam D Spatial dimensionality of the *target* operator that will use this kernel.
 *           The class itself stores a 1D kernel (inherits from @ref GaussExp\<1\>), but
 *           uses @p D to choose a normalization consistent with a D-fold separable tensor
 *           product (see notes below).
 *
 * @details
 * The constructor creates a single 1D Gaussian
 * \f[
 *   g(x) \;=\; c \, e^{-\alpha x^2},\qquad
 *   \alpha \equiv \frac{1}{\varepsilon},
 * \f]
 * then analytically differentiates it once in @c x to obtain a polynomial–Gaussian
 * (a @ref GaussPoly) and appends that single term to this expansion.
 *
 * ### Normalization and separability
 * - The coefficient is chosen as
 *   \f[
 *     c \;=\; \Big(\tfrac{\alpha}{\pi}\Big)^{D/2},
 *   \f]
 *   which corresponds to the *D-dimensional* unit-charge normalization of the isotropic
 *   Gaussian \f$ c \exp(-\alpha \lvert \mathbf r \rvert^2) \f$.
 * - When this 1D kernel is lifted to D dimensions as a separable product,
 *   MRCPP’s convolution machinery (@ref ConvolutionOperator) rescales each 1D factor
 *   by the D-th root of the magnitude of its coefficient so that the tensor product
 *   has the intended overall normalization. In effect, each axis receives
 *   \f$ (\alpha/\pi)^{1/2} \f$ and the product recovers \f$ (\alpha/\pi)^{D/2} \f$.
 *
 * ### Width control
 * The user-provided @p epsilon controls the width via \f$ \alpha = 1/\varepsilon \f$:
 * - Small \f$ \varepsilon \Rightarrow \alpha \gg 1 \Rightarrow \f$ very narrow kernel,
 *   closer to a distributional derivative, but harder to resolve numerically.
 * - Large \f$ \varepsilon \Rightarrow \alpha \ll 1 \Rightarrow \f$ broad kernel,
 *   smoother but less localized derivative approximation.
 *
 * ### Usage
 * Typically constructed internally by derivative-style convolution operators
 * (e.g., @ref DerivativeConvolution) and not used directly. If used directly,
 * pass it to a @ref ConvolutionOperator builder which will project, lift, and
 * cache the corresponding multiwavelet operator blocks.
 */
template <int D> class DerivativeKernel final : public GaussExp<1> {
public:
    /**
     * @brief Construct a rank-1 1D derivative-of-Gaussian kernel.
     *
     * @param epsilon Width control parameter; the Gaussian exponent is set to
     *                \f$ \alpha = 1/\varepsilon \f$.
     *
     * @post The expansion contains a single @ref GaussPoly term equal to
     *       \f$ \frac{d}{dx}\big[c \exp(-\alpha x^2)\big] \f$ with
     *       \f$ c = (\alpha/\pi)^{D/2} \f$.
     */
    DerivativeKernel(double epsilon)
            : GaussExp<1>() {
        // Exponent (narrowness): alpha = 1 / epsilon
        double alpha = 1.0 / epsilon;

        // D-dimensional normalization chosen up-front.
        // ConvolutionOperator later redistributes this across dimensions (D-th root per axis).
        double coef = std::pow(alpha / mrcpp::pi, D / 2.0);

        // Start from a pure 1D Gaussian g(x) = coef * exp(-alpha x^2)
        GaussFunc<1> g(alpha, coef);

        // Differentiate analytically to obtain a polynomial–Gaussian (DoG) and store it
        GaussPoly<1> dg = g.differentiate(0);

        // Single-term expansion: { dg }
        this->append(dg);
    }
};

} // namespace mrcpp
