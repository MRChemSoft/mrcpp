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

#include "ScalingBasis.h"

namespace mrcpp {

/**
 * @class InterpolatingBasis
 * @brief Interpolating multiwavelet scaling basis as defined by Alpert et al., J. Comput. Phys. 182, 149 (2002)
 *
 * @details
 * Represents the interpolatory scaling functions for the multiwavelet framework. These functions satisfy the
 * cardinal property \f$ \varphi_j(x_i) = \delta_{ij} \f$ at the Gauss quadrature nodes \f$ \{x_i\} \f$,
 * which makes projection and evaluation particularly efficient. The coefficient-to-value and
 * value-to-coefficient maps are diagonal. The constructor calls ScalingBasis(k, Interpol) to tag the family
 * and allocate storage, then fills #funcs, #quadVals, #cvMap, and #vcMap via initScalingBasis(),
 * calcQuadratureValues(), and calcCVMaps().
 *
 * @see LegendreBasis for the orthonormal (non-cardinal) alternative
 */

class InterpolatingBasis final : public ScalingBasis {
public:
    /**
     * @brief Construct an interpolating scaling basis of polynomial order @p k
     * @param k Polynomial order (typical range \f$ 1 < k < 40 \f$)
     *
     * @details Tags the basis as interpolating-family, builds the cardinal polynomials satisfying
     * \f$ \varphi_j(x_i) = \delta_{ij} \f$, sets #quadVals to the identity matrix, and assembles
     * the diagonal #cvMap and #vcMap from quadrature weights
     */
    InterpolatingBasis(int k)
            : ScalingBasis(k, Interpol) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }

private:
    /**
     * @brief Build and store the interpolating scaling functions up to degree @p k
     *
     * @details Fills #funcs with
     * \f[
     *   \varphi_j(x) = \sqrt{w_j} \sum_{m=0}^{k} \phi_m(x_j)\,\phi_m(x), \quad x \in (0,1), \quad j = 0,\ldots,k,
     * \f]
     * where \f$ \phi_m \f$ are the Legendre scaling functions and \f$ w_j \f$ are the Gauss quadrature weights.
     * The resulting functions satisfy the cardinal property \f$ \varphi_j(x_i) = \delta_{ij} \f$ at the
     * Gauss nodes \f$ \{x_i\} \f$.
     */
    void initScalingBasis();

    /**
     * @brief Set #quadVals to the identity matrix
     *
     * @details Exploits the cardinal property: evaluating \f$ \varphi_k \f$ at node \f$ x_j \f$ gives
     * \f$ \delta_{kj} \f$, so the evaluation matrix is the \f$ q \times q \f$ identity
     */
    void calcQuadratureValues();

    /**
     * @brief Assemble the diagonal coefficient-to-value and value-to-coefficient maps
     *
     * @details Sets \f$ \text{cvMap}(k,k) = 1/\sqrt{w_k} \f$ and \f$ \text{vcMap}(k,k) = \sqrt{w_k} \f$,
     * which are exact inverses of each other under the chosen normalization
     */
    void calcCVMaps();
};

} // namespace mrcpp