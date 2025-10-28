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
 * @file IdentityKernel.h
 * @brief Gaussian surrogate of the Dirac delta kernel for use in separable
 *        convolution operators.
 *
 * @details
 * This header declares @ref mrcpp::IdentityKernel, a convenience wrapper that
 * builds a one-term @ref mrcpp::GaussExp "Gaussian expansion" approximating the
 * identity (Dirac delta) kernel in \f$ \mathbb{R}^D \f$:
 * \f[
 *   \delta(x) \;\approx\; \alpha \, e^{-\beta x^2},
 * \f]
 * with parameters chosen from a requested *narrowness* (precision) \f$ \varepsilon \f$.
 * Concretely,
 * \f[
 *   \beta = \sqrt{\tfrac{1}{\varepsilon}}, \qquad
 *   \alpha = \left( \frac{\beta}{\pi} \right)^{D/2}.
 * \f]
 *
 * The resulting object is a rank-1 @ref GaussExp<1> suitable for constructing
 * separable, bandwidth-limited identity-like convolution operators; see
 * @ref mrcpp::IdentityConvolution.
 *
 * @see IdentityConvolution, ConvolutionOperator, GaussExp, GaussFunc
 */

#pragma once

#include "functions/GaussExp.h"
#include "functions/GaussFunc.h"

namespace mrcpp {

/**
 * @class IdentityKernel
 * @ingroup kernels
 * @brief Single-term Gaussian expansion approximating the Dirac delta in \f$ \mathbb{R}^D \f$.
 *
 * @tparam D Spatial dimension for the normalization of the Gaussian surrogate.
 *
 * @details
 * Constructs a one-dimensional Gaussian \f$ \alpha e^{-\beta x^2} \f$ with
 * \f$ \beta=\sqrt{1/\varepsilon} \f$ and
 * \f$ \alpha=(\beta/\pi)^{D/2} \f$,
 * then appends it to the underlying @ref GaussExp container. The parameter
 * \p epsilon controls the narrowness of the surrogate: smaller values yield
 * narrower Gaussians (closer to a true delta) but demand more resolution.
 */
template <int D> class IdentityKernel final : public GaussExp<1> {
public:
    /**
     * @brief Build a delta-like Gaussian kernel from a target narrowness \p epsilon.
     * @param epsilon Positive parameter controlling the kernel width; smaller ⇒ narrower.
     */
    IdentityKernel(double epsilon)
            : GaussExp<1>() {
        double expo = std::sqrt(1.0 / epsilon);                 // β
        double coef = std::pow(expo / mrcpp::pi, D / 2.0);      // α = (β/π)^{D/2}
        GaussFunc<1> gFunc(expo, coef);
        this->append(gFunc);
    }
};

} // namespace mrcpp