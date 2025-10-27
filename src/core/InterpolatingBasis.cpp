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
 * Implementation of the interpolating multiwavelet *scaling* basis setup.
 *
 * The goal of this class is to build a set of polynomials {I_k} such that
 * they are *interpolatory* with respect to a chosen Gaussian quadrature:
 *   • I_k evaluated at quadrature nodes (roots) forms an identity matrix.
 *   • The basis is normalized w.r.t. the quadrature weights.
 *
 * Key ingredients:
 *   - QuadratureCache supplies (roots, weights) for a given quadrature order q.
 *   - LegendrePoly(k, 2.0, 1.0) provides a scaled/shifted Legendre polynomial
 *     of degree k (the exact affine scaling is handled by LegendrePoly).
 *   - Each interpolating basis polynomial I_k is assembled as a linear
 *     combination of Legendre polynomials, then scaled by sqrt(weight_k) so
 *     that the quadrature-induced inner-product is normalized.
 *
 * Data members touched here (belonging to InterpolatingBasis):
 *   - funcs    : vector of Polynomial objects storing the scaling basis {I_k}.
 *   - quadVals : matrix of values of basis at quadrature nodes (made identity).
 *   - cvMap    : diagonal map from coefficient-space → value-space at nodes.
 *   - vcMap    : diagonal map from value-space at nodes → coefficient-space.
 */

/*
 *
 *
 *  \date June 2, 2010
 *  \author Stig Rune Jensen \n
 *          CTCC, University of Tromsø
 *
 */

#include "InterpolatingBasis.h"

#include <cmath>

#include "QuadratureCache.h"
#include "functions/LegendrePoly.h"

using namespace Eigen;

namespace mrcpp {

/**
 * @brief Build the set of interpolating scaling polynomials {I_k}.
 *
 * @details Fills
 * std::vector<Polynomial> \b funcs
 * declared in the base class
 * @ref ScalingBasis
 * with the interpolating scaling functions
 * \f[
 *    \varphi_j(x)
 *    =
 *    \sqrt{ w_j } \sum_{m = 0}^k \phi_m(x_j) \phi_m(x)
 *    , \quad
 *    x \in (0, 1)
 *    , \quad
 *    j = 0, \ldots, k
 *    ,
 * \f]
 * where \f$ \phi_m \f$ are the Legendre scaling functions.
 * Here \f$ k \f$ is \b order declared in the base class.
 * 
 * @note These interpolating scaling functions are defined on the unit interval \f$ (0, 1) \f$.

 * Procedure (for quadrature order q and scaling order s):
 *  1) Fetch Gaussian quadrature nodes (roots) and weights (wgts) of order q.
 *  2) Precompute Legendre polynomials L_0, L_1, …, L_{q-1} (scaled/shifted
 *     variant via LegendrePoly(k, 2.0, 1.0)).
 *  3) For each quadrature node k:
 *       a) Start from a copy of L_s (highest degree used for stability).
 *       b) Scale it so that I_k(roots[k]) accumulates the desired unit
 *          contribution. The factor (2*i+1) is the standard Legendre
 *          normalization multiplier that appears in expansions / projections.
 *       c) Accumulate lower-degree Legendre polynomials down to degree 0,
 *          with coefficients proportional to L_i(roots[k]) * (2*i+1).
 *       d) Finally, scale I_k by sqrt(wgts[k]) to make the quadrature-based
 *          normalization diagonal and simple (see calcCVMaps()).
 *  4) Store I_k into this->funcs.
 *
 * Remarks:
 *  - The outer loop is over nodes k, producing one cardinal/interpolatory
 *    polynomial per node.
 *  - The inner loop goes from high to low degree (q-2 … 0). The comment in
 *    the code notes that adding higher-order polys into lower-order ones is
 *    numerically undesirable, hence the chosen order of accumulation.
 */
void InterpolatingBasis::initScalingBasis() {
    int qOrder = getQuadratureOrder();   // number of quadrature points (q)
    int sOrder = getScalingOrder();      // polynomial "scaling order" (s)

    // Obtain quadrature nodes and weights of order q.
    getQuadratureCache(qc);
    const VectorXd roots = qc.getRoots(qOrder);   // size q
    const VectorXd wgts  = qc.getWeights(qOrder); // size q

    // Precompute Legendre polynomials L_k (scaled/shifted variant).
    std::vector<LegendrePoly> L_k;
    for (int k = 0; k < qOrder; k++) { L_k.push_back(LegendrePoly(k, 2.0, 1.0)); }

    // Build one interpolating polynomial I_k for each quadrature node k.
    for (int k = 0; k < qOrder; k++) {
        // Start from a copy of L_s. The comment explains the loop-order choice:
        // We avoid "adding higher-order into lower-order"; begin at top degree.
        Polynomial I_k(L_k[sOrder]);

        // Seed I_k with the value of L_s at the k-th node times (2s+1).
        // This sets up the leading contribution at node k.
        I_k *= L_k[sOrder].evalf(roots(k)) * (2.0 * sOrder + 1);

        // Accumulate lower degrees i = q-2 down to 0:
        // Each step adds val * L_i, where val depends on L_i evaluated at
        // the current node and the usual (2i+1) normalization factor.
        for (int i = qOrder - 2; i >= 0; i--) {
            double val = L_k[i].evalf(roots(k)) * (2.0 * i + 1);
            I_k.addInPlace(val, L_k[i]);
        }

        // Normalize with the square root of the quadrature weight at node k,
        // so that later the coefficient↔value maps are simple diagonal scalings.
        I_k *= std::sqrt(wgts[k]);

        // Save the constructed interpolatory scaling function for node k.
        this->funcs.push_back(I_k);
    }
}

/**
 * @brief Fill the matrix of basis values at quadrature nodes.
 *
 * For an *interpolating* basis, evaluating basis polynomial I_k at node k'
 * yields δ_{k,k'}. Therefore, the quadrature value matrix is just the identity.
 *
 * Implementation detail:
 *  - Only the diagonal entries are set to 1; all other entries remain 0
 *    (matrix presumed zero-initialized elsewhere).
 */
void InterpolatingBasis::calcQuadratureValues() {
    int q_order = getQuadratureOrder();
    for (int k = 0; k < q_order; k++) { this->quadVals(k, k) = 1.0; }
}

/**
 * @brief Build coefficient↔value diagonal maps using quadrature weights.
 *
 * The maps relate coefficient vectors in the interpolatory basis to vectors
 * of point-values at quadrature nodes, under the quadrature-induced inner
 * product:
 *
 *   - cvMap: coefficient → value map at nodes (scales by sqrt(1/w_k))
 *   - vcMap: value → coefficient map at nodes (scales by sqrt(w_k))
 *
 * With the construction in initScalingBasis(), these maps are diagonal and
 * inverse of each other.
 */
void InterpolatingBasis::calcCVMaps() {
    int q_order = getQuadratureOrder();
    getQuadratureCache(qc);
    const VectorXd &wgts = qc.getWeights(q_order);

    for (int k = 0; k < q_order; k++) {
        this->cvMap(k, k) = std::sqrt(1.0 / wgts(k)); // coeff → values
        this->vcMap(k, k) = std::sqrt(wgts(k));       // values → coeff
    }
}

} // namespace mrcpp
