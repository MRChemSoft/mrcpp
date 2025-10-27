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

/** @class LegendreBasis
 *
 * @brief Legendre scaling functions as defined by Alpert,
 * SIAM J Math Anal 24 (1), 246 (1993).
 *
 * High-level overview
 * -------------------
 * LegendreBasis represents the *Legendre scaling functions* used as a scaling
 * space in the multiwavelet framework. In contrast to an *interpolating* basis,
 * here the basis functions are (shifted/scaled) Legendre polynomials with
 * exact L² normalization. This choice leads to dense coefficient↔value maps
 * (built from evaluations at quadrature nodes), but offers orthogonality and
 * well-understood approximation properties.
 *
 * Relationship to the class hierarchy
 * -----------------------------------
 * - Inherits from @ref ScalingBasis, which provides:
 *     • storage for basis polynomials (e.g. `funcs`),
 *     • quadrature order and data,
 *     • matrices for basis evaluated at quadrature nodes (`quadVals`),
 *     • conversion maps between coefficient and nodal value spaces
 *       (`cvMap` and `vcMap`).
 *
 * What the constructor does
 * -------------------------
 * The constructor takes the polynomial order `k` (with typical bounds 1 < k < 40)
 * and:
 *   1) calls the base `ScalingBasis(k, Legendre)` to set the family/tag,
 *   2) `initScalingBasis()` to build the list of normalized Legendre polynomials
 *      up to degree `k`,
 *   3) `calcQuadratureValues()` to evaluate the basis at quadrature nodes,
 *   4) `calcCVMaps()` to assemble value→coefficient (`vcMap`) using quadrature
 *      weights and then compute coefficient→value (`cvMap`) as its inverse.
 *
 * Notes
 * -----
 * - The actual construction details are implemented in the corresponding .cpp:
 *     • `initScalingBasis()` multiplies P_k by √(2k+1) for exact normalization.
 *     • `calcQuadratureValues()` fills `quadVals(i,k) = P_k(x_i)`.
 *     • `calcCVMaps()` sets `vcMap(i,k) = P_k(x_i) * w_i` and inverts it.
 */

class LegendreBasis final : public ScalingBasis {
public:
    /** @returns New LegendreBasis object
     * @param[in] k: Polynomial order of basis, `1 < k < 40`
     *
     * Construction sequence:
     *  - `ScalingBasis(k, Legendre)` tags this as a Legendre-family scaling basis.
     *  - `initScalingBasis()` builds normalized Legendre polynomials {P_0..P_k}.
     *  - `calcQuadratureValues()` evaluates the basis at Gaussian nodes.
     *  - `calcCVMaps()` creates value↔coefficient maps using quadrature weights.
     */
    LegendreBasis(int k)
            : ScalingBasis(k, Legendre) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }

private:
    /** @brief Build and store the normalized Legendre polynomials up to degree k. */
    void initScalingBasis();
    /** @brief Fill the matrix of basis values at quadrature nodes. */
    void calcQuadratureValues();
    /** @brief Assemble value→coefficient map and its inverse (coeff→value). */
    void calcCVMaps();
};

} // namespace mrcpp