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

/** @class InterpolatingBasis
 *
 * @brief Interpolating scaling functions as defined by Alpert etal,
 * J Comp Phys 182, 149-190 (2002).
 *
 * High-level overview
 * -------------------
 * InterpolatingBasis represents the *interpolatory scaling functions* used in
 * the multiwavelet framework. These functions are constructed so that:
 *   • they interpolate at Gaussian quadrature nodes (cardinal property),
 *   • the quadrature-induced inner product is simple/diagonal,
 *   • they form the scaling space for the chosen polynomial order.
 *
 * Relationship to the hierarchy:
 *   - Inherits from ScalingBasis, which provides common functionality for
 *     scaling-function families (orders, quadrature data, storage for basis
 *     polynomials, value/coefficient maps, etc.).
 *   - The constructor finalizes initialization by calling three private
 *     helpers:
 *        1) initScalingBasis()     — build the interpolating polynomials,
 *        2) calcQuadratureValues() — fill values at quadrature nodes,
 *        3) calcCVMaps()           — build coefficient↔value diagonal maps.
 *
 * Mathematical context (very short):
 *   - Follows the construction in Alpert (2002) for interpolatory multiwavelets,
 *     where basis functions {I_k} satisfy I_k(x_j) = δ_{k,j} at quadrature nodes
 *     {x_j}. This makes projection/evaluation particularly efficient.
 */

class InterpolatingBasis final : public ScalingBasis {
public:
    /** @returns New InterpolatingBasis object
     * @param[in] k: Polynomial order of basis, `1 < k < 40`
     *
     * What happens in the constructor:
     *  - Calls the ScalingBasis base constructor with (k, Interpol), which
     *    sets the family/type to “Interpolating”.
     *  - initScalingBasis(): constructs the set of interpolating polynomials
     *    (stored in the base's internal container, typically `funcs`).
     *  - calcQuadratureValues(): sets the basis evaluation matrix at nodes to
     *    the identity (cardinality property).
     *  - calcCVMaps(): builds diagonal conversion maps between coefficient
     *    vectors and values at quadrature nodes using the quadrature weights.
     *
     * Precondition:
     *  - k must be within the supported range of the library (checked by the
     *    base class). Typical limits are 1 < k < 40 as noted here.
     */
    InterpolatingBasis(int k)
            : ScalingBasis(k, Interpol) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }

private:
    /**
     * @brief Construct the interpolatory scaling polynomials {I_k}.
     *
     * Implementation details (in .cpp):
     *  - Uses Gaussian quadrature roots/weights of order q.
     *  - Expands I_k in a Legendre polynomial basis and enforces I_k(x_j)=δ_{kj}.
     *  - Applies sqrt(weight) normalization so that the induced inner product
     *    is diagonal and the cv/vc maps become simple scalings.
     */
    void initScalingBasis();

    /**
     * @brief Fill the basis-at-nodes matrix.
     *
     * For an interpolating basis, evaluating the k-th basis at node j yields
     * δ_{kj}. The implementation sets the diagonal entries to 1 (identity).
     */
    void calcQuadratureValues();

    /**
     * @brief Build coefficient↔value diagonal maps using quadrature weights.
     *
     * - cvMap(k,k) = sqrt(1 / w_k)  (coefficients → values at nodes)
     * - vcMap(k,k) = sqrt(w_k)      (values at nodes → coefficients)
     *
     * These maps are exact inverses due to the chosen normalization.
     */
    void calcCVMaps();
};

} // namespace mrcpp