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
#include "eigen_disable_warnings.h"

using namespace Eigen;

namespace mrcpp {

void LegendreBasis::initScalingBasis() {
    for (int k = 0; k < getScalingOrder() + 1; k++) {
        LegendrePoly L_k(k, 2.0, 1.0);
        L_k *= std::sqrt(2.0 * k + 1.0); // exact normalization
        this->funcs.push_back(L_k);
    }
}

void LegendreBasis::calcQuadratureValues() {
    getQuadratureCache(qc);
    int q_order = getQuadratureOrder();
    const VectorXd &pts = qc.getRoots(q_order);

    for (int k = 0; k < q_order; k++) {
	const Polynomial &poly = this->getFunc(k);
	for (int i = 0; i < q_order; i++) {
	    this->quadVals(i, k) = poly.evalf(pts(i));
	}
    }
}

void LegendreBasis::calcCVMaps() {
    getQuadratureCache(qc);
    int q_order = getQuadratureOrder();
    const VectorXd &pts = qc.getRoots(q_order);
    const VectorXd &wgts = qc.getWeights(q_order);

    for (int k = 0; k < q_order; k++) {
	const Polynomial &poly = this->getFunc(k);
	for (int i = 0; i < q_order; i++) {
	    this->vcMap(i, k) = poly.evalf(pts(i)) * wgts(i);
	}
    }
    this->cvMap = this->vcMap.inverse();
}

} // namespace mrcpp
