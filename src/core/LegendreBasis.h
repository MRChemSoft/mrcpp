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

#include <Eigen/Dense>

#include "ScalingBasis.h"

namespace mrcpp {

/**
 * @class LegendreBasis
 * @brief Legendre multiwavelet scaling basis as defined by Alpert, SIAM J. Math. Anal. 24(1), 246 (1993)
 *
 * @details
 * Represents the Legendre scaling functions used as the scaling space in the multiwavelet framework.
 * The basis functions are orthonormalized shifted Legendre polynomials:
 * \f[
 *   \varphi_j(x) = \sqrt{2j+1}\,P_j(2x-1), \quad x \in (0,1), \quad j = 0,\ldots,k.
 * \f]
 * In contrast to the interpolating basis, the coefficient-to-value map (built from evaluations at
 * quadrature nodes) is dense. The constructor calls the base ScalingBasis(k, Legendre) to tag the
 * family and size the matrices, then populates #funcs, #quadVals, #cvMap, and #vcMap via the three
 * private helpers initScalingBasis(), calcQuadratureValues(), and calcCVMaps().
 *
 * @see InterpolatingBasis for the cardinal (interpolatory) alternative
 */

class LegendreBasis final : public ScalingBasis {
public:
    /**
     * @brief Construct a Legendre scaling basis of polynomial order @p k
     * @param k Polynomial order (typical range \f$ 1 < k < 40 \f$)
     *
     * @details Tags the basis as Legendre-family, builds the \f$ q = k+1 \f$ normalized Legendre
     * polynomials, evaluates them at Gauss nodes to fill #quadVals, and assembles #vcMap and #cvMap
     */
    LegendreBasis(int k)
            : ScalingBasis(k, Legendre) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }

private:
    /**
     * @brief Build and store the normalized Legendre scaling functions up to degree k
     *
     * @details Fills the base-class container @c funcs with
     * \f[
     *   \phi_j(x) = \sqrt{2j+1}\,P_j(2x-1),
     *   \quad x \in (0,1),\quad j = 0,\ldots,k,
     * \f]
     * where \f$ P_j \f$ are standard Legendre polynomials and \f$ k \f$ is the
     * polynomial order declared in @ref ScalingBasis
     *
     * @note The Legendre scaling functions are defined on the unit interval \f$ (0,1) \f$
     */
    void initScalingBasis();
    /** @brief Fill the matrix of basis values at quadrature nodes */
    void calcQuadratureValues();
    /** @brief Assemble value→coefficient map and its inverse (coeff→value) */
    void calcCVMaps();
};

} // namespace mrcpp