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
#include "LegendrePoly.h"

using namespace std;
using namespace Eigen;

void InterpolatingBasis::initScalingBasis() {
    int qOrder = getQuadratureOrder();
    int sOrder = getScalingOrder();

    getQuadratureCache(qc);
    const VectorXd roots = qc.getRoots(qOrder);
    const VectorXd wgts = qc.getWeights(qOrder);

    vector<LegendrePoly> L_k;
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
        I_k *= sqrt(wgts[k]);
        this->funcs.push_back(I_k);
    }
}

//void InterpolatingBasis::preEvaluate() {
//	getQuadratureCache(qc);
//	const VectorXd &wgts = qc.getWeights(quadratureOrder);

//	int npts = scalingOrder + 1;
//	preVals = MatrixXd::Zero(quadratureOrder, npts);

//	assert(quadratureOrder == npts);

//	for (int k = 0; k < npts; k++) {
//		preVals(k, k) = sqrt(wgts(k));
//	}
//}

///****** WARNING! Ugliness ahead!!! ********************************/

//VectorXd InterpolatingBasis::calcScalingCoefs(int d,
//	const SeparableFunction<1> &func, int n, int l) const {
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
//		scs.col(d)(i) = func.evalf(r, d) * preVals(i, i);
//	}
//	scs *= sqrt(1.0 / sfac);
//	return scs;
//}

//VectorXd InterpolatingBasis::calcScalingCoefs(int d,
//	const SeparableFunction<2> &func, int n, int l) const {
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
//		scs.col(d)(i) = func.evalf(r, d) * preVals(i, i);
//	}
//	scs *= sqrt(1.0 / sfac);
//	return scs;
//}

//VectorXd InterpolatingBasis::calcScalingCoefs(int d,
//	const SeparableFunction<3> &func, int n, int l) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);
//	VectorXd scs(scalingOrder + 1);

//	double r;
//	double sfac;
//	if (n < 0) {
//		sfac = pow(2.0, n);
//	} else {
//		sfac = 1 << n;
//	}

//	for (int i = 0; i < quadratureOrder; i++) {
//		r = (pts(i) + l) / sfac;
//		scs.col(d)(i) = func.evalf(r, d) * preVals(i, i);
//	}
//	scs *= sqrt(1.0 / sfac);
//	return scs;
//}

//void InterpolatingBasis::calcScalingCoefs(const SeparableFunction<1> &func,
//		int n, const int *l, MatrixXd &scs) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);
//	int dim = scs.cols();

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
//			scs.col(d)(i) = func.evalf(r, d) * preVals(i,i);
//		}
//	}
//	scs *= sqrt(1.0 / sfac);
//}

//void InterpolatingBasis::calcScalingCoefs(const SeparableFunction<2> &func,
//		int n, const int *l, MatrixXd &scs) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);
//	int dim = scs.cols();

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
//			scs.col(d)(i) = func.evalf(r, d) * preVals(i,i);
//		}
//	}
//	scs *= sqrt(1.0 / sfac);
//}

//void InterpolatingBasis::calcScalingCoefs(const SeparableFunction<3> &func,
//		int n, const int *l, MatrixXd &scs) const {
//	getQuadratureCache(qc);
//	const VectorXd &pts = qc.getRoots(quadratureOrder);
//	int dim = scs.cols();

//	double r;
//	double sfac = 1.0;
//	if (n < 0) {
//		sfac = pow(2.0, n);
//	} else if (n > 0) {
//		sfac = 1 << n;
//	}

//	for (int d = 0; d < dim; d++) {
//		for (int i = 0; i < quadratureOrder; i++) {
//			r = (pts(i) + l[d]) / sfac;
//			scs.col(d)(i) = func.evalf(r, d) * preVals(i,i);
//		}
//	}
//	if (n != 0) {
//		scs *= sqrt(1.0 / sfac);
//	}
//}

