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
 * @brief Initialise the interpolating scaling basis and populate `funcs`.
 *
 * @details
 * This builds the set of interpolating (cardinal) scaling polynomials
 * \f[
 *   \varphi_j(x)
 *   \;=\;
 *   \sqrt{w_j}\;\sum_{m=0}^{k} \phi_m(x_j)\,\phi_m(x),
 *   \qquad x\in(0,1),\; j=0,\ldots,k,
 * \f]
 * where \f$\{\phi_m\}\f$ are the Legendre scaling functions, \f$\{x_j,w_j\}\f$
 * are the Gauss–Legendre nodes and weights of order \f$q=k+1\f$, and
 * \f$k=\texttt{order}\f$ from the base class. The resulting functions are
 * stored in `std::vector<Polynomial> funcs` inherited from @ref ScalingBasis .
 *
 * Procedural outline (q := quadrature order, s := scaling order = k):
 *  1) Fetch q-point Gaussian quadrature nodes `roots` and weights `wgts`.
 *  2) Precompute Legendre polynomials \f$L_0,\ldots,L_{q-1}\f$ on (0,1) using
 *     `LegendrePoly(i, 2.0, 1.0)` (scaled/shifted variant).
 *  3) For each quadrature node j:
 *       a) Start from a copy of \f$L_s\f$ (highest degree used for stability).
 *       b) Scale so that \f$I_j(x_j)\f$ accumulates the desired unit
 *          contribution; the factor \f$(2i+1)\f$ is the standard Legendre
 *          normalization used in projections/expansions.
 *       c) Accumulate lower-degree Legendre polynomials down to degree 0 with
 *          coefficients \f$L_i(x_j)\,(2i+1)\f$.
 *       d) Multiply the resulting polynomial by \f$\sqrt{w_j}\f$ so that the
 *          quadrature-based normalization becomes diagonal (see `calcCVMaps()`).
 *  4) Push the constructed \f$I_j\f$ into `this->funcs`.
 *
 * Remarks:
 *  - One interpolatory polynomial \f$I_j\f$ is produced per node \f$x_j\f$.
 *  - The accumulation loop goes from high to low degree (\f$q-2,\ldots,0\f$)
 *    to avoid adding higher-order polynomials into lower-order ones, which is
 *    numerically less stable.
 *  - These interpolating scaling functions are defined on the unit interval
 *    \f$(0,1)\f$.
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
 * For an interpolating basis the cardinal property holds:
 * I_k(x_{k'}) = δ_{k,k'}.
 * Therefore `quadVals` is simply the identity matrix of size q×q,
 * where q = getQuadratureOrder().
 *
 * Implementation detail:
 * - Only the diagonal entries are set to 1.0; all off-diagonals remain 0.
 *   The matrix is assumed to have been zero-initialized.
 */
void InterpolatingBasis::calcQuadratureValues() {
    int q_order = getQuadratureOrder();
    for (int k = 0; k < q_order; k++) { this->quadVals(k, k) = 1.0; }
}
/**
 * @brief Build coefficient↔value diagonal maps using quadrature weights.
 *
 * These maps relate coefficient vectors in the interpolatory basis to vectors
 * of point-values at quadrature nodes under the quadrature-induced inner product.
 *
 *  - cvMap (coeff → values at nodes):   cvMap(k,k) = sqrt(1 / w_k)
 *  - vcMap (values at nodes → coeff):   vcMap(k,k) = sqrt(w_k)
 *
 * Because the basis is cardinal (I_k(x_{k'}) = δ_{k,k'}), both maps are diagonal
 * and are mutual inverses. Here {w_k} are the Gauss–Legendre weights of order q.
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