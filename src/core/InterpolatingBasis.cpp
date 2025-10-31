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

#include "InterpolatingBasis.h"

#include <cmath>

#include "QuadratureCache.h"
#include "functions/LegendrePoly.h"

using namespace Eigen;

namespace mrcpp {

void InterpolatingBasis::initScalingBasis() {
    int qOrder = getQuadratureOrder();
    int sOrder = getScalingOrder();

    getQuadratureCache(qc);
    const VectorXd roots = qc.getRoots(qOrder);
    const VectorXd wgts  = qc.getWeights(qOrder);

    std::vector<LegendrePoly> L_k;
    for (int k = 0; k < qOrder; k++) { L_k.push_back(LegendrePoly(k, 2.0, 1.0)); }

    for (int k = 0; k < qOrder; k++) {
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

void InterpolatingBasis::calcQuadratureValues() {
    int q_order = getQuadratureOrder();
    for (int k = 0; k < q_order; k++) { this->quadVals(k, k) = 1.0; }
}

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