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

namespace mrcpp {

/**
 * @class HeatKernel
 * @ingroup functions
 *
 * @brief Single-term Gaussian expansion that represents the \(D\)-dimensional
 *        heat kernel \(K_t(\mathbf x)\) at diffusion time \(t>0\).
 *
 * @tparam D Spatial dimension of the kernel to be modeled (1, 2, or 3).
 *
 * @details
 * The continuous heat kernel in \(\mathbb R^D\) is
 * \f[
 *   K_t(\mathbf x)
 *   \;=\;
 *   \frac{1}{(4\pi t)^{D/2}}
 *   \exp\!\left(-\frac{\lVert \mathbf x\rVert^2}{4t}\right),
 *   \qquad t>0.
 * \f]
 *
 * In MRCPP, separable operators are commonly assembled from 1D building blocks.
 * This class therefore inherits from @ref GaussExp "GaussExp<1>" and stores a
 * *single* 1D Gaussian term whose exponent and coefficient are chosen so that,
 * when used inside separable constructions (e.g., @ref ConvolutionOperator),
 * the resulting operator corresponds to the \(D\)-dimensional heat kernel.
 *
 * Concretely, with
 * \f$ \beta = \frac{1}{4t} \f$ and \f$ \alpha = \big(\beta/\pi\big)^{D/2} \f$,
 * we append the 1D Gaussian
 * \f[
 *   g(x) = \alpha\, e^{-\beta x^2},
 * \f]
 * and the higher-dimensional operator logic (tensor products across coordinates)
 * recovers the isotropic \(D\)-dimensional kernel.
 *
 * ### Notes
 * - The constructor does not enforce \(t>0\) at runtime; pass a strictly positive
 *   value to avoid nonsensical parameters.
 * - The class is intentionally minimal: it only sets up the Gaussian parameters
 *   and leaves projection/assembly to the caller (e.g., convolution operators).
 *
 * ### Example
 * @code
 * MultiResolutionAnalysis<3> mra(box, basis);
 * HeatKernel<3> Kt(0.05); // 3D heat kernel at t = 0.05
 * // Use Kt as a kernel for a ConvolutionOperator<3>, etc.
 * @endcode
 */
template <int D> class HeatKernel final : public GaussExp<1> {
public:
    /**
     * @brief Construct a heat kernel at diffusion time @p t.
     *
     * @param t Diffusion time (\f$ t>0 \f$). Smaller values yield narrower
     *          Gaussians (more localized kernels).
     *
     * @details
     * Sets the Gaussian exponent to \f$ \beta = \frac{1}{4t} \f$ and the
     * coefficient to \f$ \alpha = \big(\beta/\pi\big)^{D/2} \f$, then appends a
     * single @ref GaussFunc "GaussFunc<1>" to this @ref GaussExp "GaussExp<1>".
     */
    HeatKernel(double t)
            : GaussExp<1>() {
        // Exponent β = 1/(4t)
        double expo = 0.25 / t;

        // Amplitude α = (β/π)^{D/2} so that the separable product matches (4πt)^{-D/2}
        double coef = std::pow(expo / mrcpp::pi, D / 2.0);

        // Build the 1D Gaussian term and register it in the expansion
        GaussFunc<1> gFunc(expo, coef);
        this->append(gFunc);
    }
};

} // namespace mrcpp