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
 * @file HelmholtzKernel.h
 * @brief Declaration of a Gaussian expansion approximating the 3D Helmholtz (screened Coulomb) kernel.
 *
 * @details
 * This header declares @c HelmholtzKernel, a convenience wrapper that builds a
 * separable Gaussian expansion
 * \f[
 *   K_\mu(r) \;\approx\; \sum_{m=1}^{M} \beta_m\, e^{-\alpha_m r^2},
 * \f]
 * that approximates the radial 3D Helmholtz/Yukawa kernel on a finite interval
 * \f$[r_\min,r_\max]\f$ with a target relative accuracy \f$\varepsilon\f$.
 * The class derives from @ref mrcpp::GaussExp "GaussExp<1>" and therefore can be
 * used anywhere a one–dimensional Gaussian expansion is expected (e.g. to form
 * convolution operators).
 */

#pragma once

#include "functions/GaussExp.h"

namespace mrcpp {

/**
 * @class HelmholtzKernel
 * @brief Gaussian expansion of the 3D Helmholtz (screened Coulomb / Yukawa) kernel.
 *
 * @details
 * Constructs a 1D Gaussian expansion (in the radial variable) by sampling an
 * integral representation of the Helmholtz kernel in a logarithmic parameter and
 * applying a trapezoidal quadrature. The resulting set of Gaussian terms
 * \f$\{\alpha_m,\beta_m\}\f$ is rescaled to the requested physical interval
 * \f$[r_\min,r_\max]\f$.
 *
 * Typical usage:
 * @code
 *   double mu = 1.0;          // screening parameter
 *   double eps = 1e-8;        // target relative accuracy
 *   double rmin = 1e-3, rmax = 10.0;
 *   mrcpp::HelmholtzKernel kernel(mu, eps, rmin, rmax);
 *   // 'kernel' is a GaussExp<1> and can be used to build convolution operators
 * @endcode
 *
 * @note The actual separation rank @f$M@f$ depends on @p epsilon and the interval
 *       size. Extremely tight tolerances or very wide intervals may require a rank
 *       larger than the internal limit (see @c MaxSepRank in the implementation).
 */
class HelmholtzKernel final : public GaussExp<1> {
public:
    /**
     * @brief Build a Gaussian expansion of the Helmholtz kernel on \f$[r_\min,r_\max]\f$.
     *
     * @param mu       Screening parameter \f$\mu > 0\f$ (Yukawa wavenumber).
     * @param epsilon  Target relative accuracy \f$0 < \varepsilon < 1\f$ for the expansion.
     * @param r_min    Lower radius bound (must satisfy \f$0 < r_\min < r_\max\f$).
     * @param r_max    Upper radius bound.
     *
     * @details
     * The constructor fills the underlying @ref GaussExp "GaussExp<1>" with
     * \f$M\f$ Gaussian terms determined by a trapezoidal discretization in a
     * logarithmic variable. Endpoints are weighted with half-quadrature weights.
     */
    HelmholtzKernel(double mu, double epsilon, double r_min, double r_max);
};

} // namespace mrcpp
