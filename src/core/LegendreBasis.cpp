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

using namespace std;
using namespace Eigen;

namespace mrcpp {

void LegendreBasis::initScalingBasis() {
    for (int k = 0; k < getScalingOrder() + 1; k++) {
        LegendrePoly L_k(k, 2.0, 1.0);
        L_k *= sqrt(2.0 * k + 1.0); // exact normalization
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

///****** WARNING! Ugliness ahead!!! ********************************/

//VectorXd LegendreBasis::calcScalingCoefs(int d,
//		const SeparableFunction<1> &func, int n, int l) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);

//	VectorXd scs(scalingOrder + 1);
//	RowVectorXd fvals(quadratureOrder);

//	double r;
//	double sfac;
//	if (n < 0) {
//		sfac = pow(2.0, n);
//	} else {
//		sfac = 1 << n;
//	}

//	for (int i = 0; i < quadratureOrder; i++) {
//		r = (pts(i) + l) / sfac;
//		fvals(i) = func.evalf(r, d);
//	}
//	scs = fvals * preVals;
//	scs *= sqrt(1.0 / sfac);
//	return scs;
//}

//VectorXd LegendreBasis::calcScalingCoefs(int d,
//		const SeparableFunction<2> &func, int n, int l) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);

//	VectorXd scs(scalingOrder + 1);
//	RowVectorXd fvals(quadratureOrder);

//	double r;
//	double sfac;
//	if (n < 0) {
//		sfac = pow(2.0, n);
//	} else {
//		sfac = 1 << n;
//	}

//	for (int i = 0; i < quadratureOrder; i++) {
//		r = (pts(i) + l) / sfac;
//		fvals(i) = func.evalf(r, d);
//	}
//	scs = fvals * preVals;
//	scs *= sqrt(1.0 / sfac);
//	return scs;
//}

//VectorXd LegendreBasis::calcScalingCoefs(int d,
//		const SeparableFunction<3> &func, int n, int l) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);

//	VectorXd scs(scalingOrder + 1);
//	RowVectorXd fvals(quadratureOrder);

//	double r;
//	double sfac;
//	if (n < 0) {
//		sfac = pow(2.0, n);
//	} else {
//		sfac = 1 << n;
//	}

//	for (int i = 0; i < quadratureOrder; i++) {
//		r = (pts(i) + l) / sfac;
//		fvals(i) = func.evalf(r, d);
//	}
//	scs = fvals * preVals;
//	scs *= sqrt(1.0 / sfac);
//	return scs;
//}

//void LegendreBasis::calcScalingCoefs(const SeparableFunction<1> &func,
//		int n, const int *l, MatrixXd &scs) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);

//	int dim = scs.cols();
//	RowVectorXd fvals(quadratureOrder);

//	double r;
//	double sfac;
//	if (n < 0) {
//		sfac = pow(2.0, n);
//	} else {
//		sfac = 1 << n;
//	}

//	for (int d = 0; d < dim; d++) {
//		for (int i = 0; i < quadratureOrder; i++) {
//			r = (pts(i) + l[d]) / sfac;
//			fvals(i) = func.evalf(r, d);
//		}
//		scs.col(d) = fvals * preVals;
//	}

//	scs *= sqrt(1.0 / sfac);
//}

//void LegendreBasis::calcScalingCoefs(const SeparableFunction<2> &func,
//		int n, const int *l, MatrixXd &scs) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);

//	int dim = scs.cols();
//	RowVectorXd fvals(quadratureOrder);

//	double r;
//	double sfac;
//	if (n < 0) {
//		sfac = pow(2.0, n);
//	} else {
//		sfac = 1 << n;
//	}

//	for (int d = 0; d < dim; d++) {
//		for (int i = 0; i < quadratureOrder; i++) {
//			r = (pts(i) + l[d]) / sfac;
//			fvals(i) = func.evalf(r, d);
//		}
//		scs.col(d) = fvals * preVals;
//	}

//	scs *= sqrt(1.0 / sfac);
//}

//void LegendreBasis::calcScalingCoefs(const SeparableFunction<3> &func,
//		int n, const int *l, MatrixXd &scs) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);

//	int dim = scs.cols();
//	RowVectorXd fvals(quadratureOrder);

//	double r;
//	double sfac;
//	if (n < 0) {
//		sfac = pow(2.0, n);
//	} else {
//		sfac = 1 << n;
//	}

//	for (int d = 0; d < dim; d++) {
//		for (int i = 0; i < quadratureOrder; i++) {
//			r = (pts(i) + l[d]) / sfac;
//			fvals(i) = func.evalf(r, d);
//		}
//		scs.col(d) = fvals * preVals;
//	}
//	scs *= sqrt(1.0 / sfac);
//}

} // namespace mrcpp
