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


/** @brief Initialise interpolating scaling basis.
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
 * 
 *
 * 
 */
void InterpolatingBasis::initScalingBasis() {
    int qOrder = getQuadratureOrder();
    int sOrder = getScalingOrder();     // sOrder = qOrder - 1

    getQuadratureCache(qc);
    const VectorXd roots = qc.getRoots(qOrder);
    const VectorXd wgts = qc.getWeights(qOrder);

    std::vector<LegendrePoly> L_k;
    for (int k = 0; k < qOrder; k++) { L_k.push_back(LegendrePoly(k, 2.0, 1.0)); }

    for (int k = 0; k < qOrder; k++) {
        // Can't add higher-order polynomials to lower-order ones, so I
        // changed the order of the loop
        Polynomial I_k(L_k[sOrder]);
        I_k *= L_k[sOrder].evalf(roots(k)) * (2.0 * sOrder + 1);

        for (int i = qOrder - 2; i >= 0; i--) {
            double val = L_k[i].evalf(roots(k)) * (2.0 * i + 1);
            I_k.addInPlace(val, L_k[i]);
        }
        I_k *= std::sqrt(wgts[k]);
        this->funcs.push_back(I_k);
    }
}


/** @brief In Progress by Evgueni...
 * 
 *
 * 
 */
void InterpolatingBasis::calcQuadratureValues() {
    int q_order = getQuadratureOrder();
    for (int k = 0; k < q_order; k++) { this->quadVals(k, k) = 1.0; }
}


/** @brief In Progress by Evgueni...
 * 
 *
 * 
 */
void InterpolatingBasis::calcCVMaps() {
    int q_order = getQuadratureOrder();
    getQuadratureCache(qc);
    const VectorXd &wgts = qc.getWeights(q_order);

    for (int k = 0; k < q_order; k++) {
        this->cvMap(k, k) = std::sqrt(1.0 / wgts(k));
        this->vcMap(k, k) = std::sqrt(wgts(k));
    }
}

} // namespace mrcpp
