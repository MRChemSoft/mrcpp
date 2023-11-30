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

#include "LegendreBasis.h"
#include "QuadratureCache.h"
#include "functions/LegendrePoly.h"

using namespace Eigen;

namespace mrcpp {


/** @brief Initialise Legendre scaling basis.
 * 
 * @details Fills
 * std::vector<Polynomial> \b funcs
 * declared in the base class
 * @ref ScalingBasis
 * with the Legendre scaling functions
 * \f[
 *    \phi_j(x)
 *    =
 *    \sqrt{ 2j + 1 } P_j(2x - 1)
 *    , \quad
 *    x \in (0, 1)
 *    , \quad
 *    j = 0, \ldots, k
 *    ,
 * \f]
 * where \f$ P_j \f$ are standard Legendre polynomials.
 * Here \f$ k \f$ is \b order declared in the base class.
 * 
 * @note These Legendre scaling functions are defined on the unit interval \f$ (0, 1) \f$.
 * 
 */
void LegendreBasis::initScalingBasis() {
    for (int k = 0; k < getScalingOrder() + 1; k++) {
        LegendrePoly L_k(k, 2.0, 1.0);
        L_k *= std::sqrt(2.0 * k + 1.0); // exact normalization
        this->funcs.push_back(L_k);
    }
}


/** @brief In Progress by Evgueni...
 * 
 *
 * 
 */
void LegendreBasis::calcQuadratureValues() {
    getQuadratureCache(qc);
    int q_order = getQuadratureOrder();
    const VectorXd &pts = qc.getRoots(q_order);

    for (int k = 0; k < q_order; k++) {
        const Polynomial &poly = this->getFunc(k);
        for (int i = 0; i < q_order; i++) { this->quadVals(i, k) = poly.evalf(pts(i)); }
    }
}


/** @brief In Progress by Evgueni...
 * 
 *
 * 
 */
void LegendreBasis::calcCVMaps() {
    getQuadratureCache(qc);
    int q_order = getQuadratureOrder();
    const VectorXd &pts = qc.getRoots(q_order);
    const VectorXd &wgts = qc.getWeights(q_order);

    for (int k = 0; k < q_order; k++) {
        const Polynomial &poly = this->getFunc(k);
        for (int i = 0; i < q_order; i++) { this->vcMap(i, k) = poly.evalf(pts(i)) * wgts(i); }
    }
    this->cvMap = this->vcMap.inverse();
}

} // namespace mrcpp
