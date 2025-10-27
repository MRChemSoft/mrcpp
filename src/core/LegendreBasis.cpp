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

/*
 * Overview
 * --------
 * Implementation of the *Legendre* scaling basis used by the multiwavelet
 * framework. In contrast to the interpolating basis, this basis consists of
 * (shifted/scaled) Legendre polynomials with exact L^2-normalization.
 *
 * Responsibilities of this file:
 *  - Build the list of scaling polynomials {P_k} up to the scaling order.
 *  - Evaluate these polynomials at Gaussian quadrature nodes to populate
 *    the quadrature value matrix (basis-at-nodes).
 *  - Construct the coefficient↔value maps using quadrature weights; here
 *    vcMap is assembled directly and cvMap is its matrix inverse.
 *
 * Notation:
 *  - LegendrePoly(k, 2.0, 1.0) represents the degree-k Legendre polynomial
 *    evaluated on an affine-mapped interval (handled by LegendrePoly).
 *  - getScalingOrder() returns the polynomial order "s".
 *  - getQuadratureOrder() returns the number of quadrature nodes "q".
 *  - funcs    : container of basis polynomials (in the base class).
 *  - quadVals : matrix of basis values at quadrature nodes (size q×(s+1)).
 *  - vcMap    : value→coefficient map built from basis values and weights.
 *  - cvMap    : inverse of vcMap (coefficient→value).
 */

/*
 *
 *
 *  \date June 2, 2010
 *  \author Stig Rune Jensen \n
 *          CTCC, University of Tromsø
 *
 */

#include "LegendreBasis.h"
#include "QuadratureCache.h"
#include "functions/LegendrePoly.h"

using namespace Eigen;

namespace mrcpp {

/**
 * @brief Initialize the Legendre scaling basis {P_k}_{k=0..s}.
 *
 * For each degree k up to the scaling order, construct a (shifted/scaled)
 * Legendre polynomial and multiply by sqrt(2k+1) to achieve exact L^2
 * normalization on the reference interval used by LegendrePoly.
 *
 * Effects:
 *  - Appends each normalized polynomial to this->funcs.
 */
void LegendreBasis::initScalingBasis() {
    for (int k = 0; k < getScalingOrder() + 1; k++) {
        LegendrePoly L_k(k, 2.0, 1.0);                 // degree-k Legendre (mapped)
        L_k *= std::sqrt(2.0 * k + 1.0);               // exact normalization factor
        this->funcs.push_back(L_k);                    // store in basis list
    }
}

/**
 * @brief Fill the matrix of basis values at Gaussian quadrature points.
 *
 * quadVals(i, k) := P_k( x_i ),  where {x_i} are the q quadrature nodes.
 *
 * Steps:
 *  1) Obtain quadrature roots (points) of order q.
 *  2) For each basis polynomial P_k, evaluate it at all points x_i and store
 *     in the corresponding column k of quadVals.
 */
void LegendreBasis::calcQuadratureValues() {
    getQuadratureCache(qc);
    int q_order = getQuadratureOrder();
    const VectorXd &pts = qc.getRoots(q_order);        // x_i, i = 0..q-1

    for (int k = 0; k < q_order; k++) {
        const Polynomial &poly = this->getFunc(k);     // P_k
        for (int i = 0; i < q_order; i++) {
            this->quadVals(i, k) = poly.evalf(pts(i)); // quadVals(i,k) = P_k(x_i)
        }
    }
}

/**
 * @brief Build the coefficient↔value maps using quadrature weights.
 *
 * For the Legendre basis, we assemble vcMap directly via:
 *   vcMap(i, k) = P_k( x_i ) * w_i,
 * where {w_i} are the quadrature weights. This corresponds to the (discrete)
 * projection of the basis onto the quadrature nodes with weighting.
 *
 * Then we compute cvMap as the matrix inverse of vcMap:
 *   cvMap = (vcMap)^{-1}.
 *
 * Interpretation:
 *  - vcMap : value→coefficient (takes nodal values and produces coefficients)
 *  - cvMap : coefficient→value (evaluates coefficients back to nodal values)
 *
 * Note:
 *  - Unlike the interpolating basis (where maps are diagonal), for the
 *    Legendre basis vcMap is dense (q×q) and we invert it numerically.
 */
void LegendreBasis::calcCVMaps() {
    getQuadratureCache(qc);
    int q_order = getQuadratureOrder();
    const VectorXd &pts  = qc.getRoots(q_order);       // x_i
    const VectorXd &wgts = qc.getWeights(q_order);     // w_i

    // Assemble vcMap(i,k) = P_k(x_i) * w_i
    for (int k = 0; k < q_order; k++) {
        const Polynomial &poly = this->getFunc(k);
        for (int i = 0; i < q_order; i++) {
            this->vcMap(i, k) = poly.evalf(pts(i)) * wgts(i);
        }
    }

    // Invert to obtain cvMap (coefficient→value).
    this->cvMap = this->vcMap.inverse();
}

} // namespace mrcpp