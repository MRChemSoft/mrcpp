/*
 *
 *
 *  \date June 2, 2010
 *  \author Stig Rune Jensen \n
 *          CTCC, University of Tromsø
 *
 */

#include "InterpolatingBasis.h"
#include "QuadratureCache.h"
#include "functions/LegendrePoly.h"

using namespace Eigen;

namespace mrcpp {

void InterpolatingBasis::initScalingBasis() {
    int qOrder = getQuadratureOrder();
    int sOrder = getScalingOrder();

    getQuadratureCache(qc);
    const VectorXd roots = qc.getRoots(qOrder);
    const VectorXd wgts = qc.getWeights(qOrder);

    std::vector<LegendrePoly> L_k;
    for (int k = 0; k < qOrder; k++) {
        L_k.push_back(LegendrePoly(k, 2.0, 1.0));
    }

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

void InterpolatingBasis::calcQuadratureValues() {
    int q_order = getQuadratureOrder();
    for (int k = 0; k < q_order; k++) {
        this->quadVals(k, k) = 1.0;
    }
}

void InterpolatingBasis::calcCVMaps() {
    int q_order = getQuadratureOrder();
    getQuadratureCache(qc);
    const VectorXd &wgts = qc.getWeights(q_order);

    for (int k = 0; k < q_order; k++) {
        this->cvMap(k, k) = std::sqrt(1.0/wgts(k));
        this->vcMap(k, k) = std::sqrt(wgts(k));
    }
}

} // namespace mrcpp
