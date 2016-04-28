/**
 * 
 *
 *  \date Jul 7, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include <iostream>

#include "GreensKernel.h"
#include "BoundingBox.h"

using namespace std;
double GreensKernel::defaultRMin = 1.0e-5;
double GreensKernel::defaultRMax = 1.75;
int GreensKernel::defaultCompMin = 0;
int GreensKernel::defaultCompMax = MAX_SEP_RANK;
double GreensKernel::defaultExpPrec = 1.0e-5;
double GreensKernel::defaultProjPrec = 1.0e-6;

/** convolution of a kernel defined as a multidim gaussian function with
 *  a function defined as multidim gaussian expansion \f$ g = k * f \f$
 */
GaussExp<1> GreensKernel::apply(Gaussian<1> &f) {
	return convolute(f);
}

GaussExp<1> GreensKernel::apply(GaussExp<1> &f) {
	int nfterms = f.size();
	int nkterms = kern.size();
	GaussExp<1> convl(nfterms * nkterms);

	int n = 0;
	for (int i = 0; i < nfterms; i++) {
		for (int k = 0; k < nkterms; k++) {
			const Gaussian<1> *tmp=&f.getFunc(i);
			GaussFunc<1> g = convolute1(k, *tmp);
			convl.setFunc(n, g);
			n++;
		}
	}
	return convl;
}

/** convolution of a one-dimensional kernel with a function defined as
 *  multidim gaussian expansion \f$ g = k * f \f$
 */
GaussFunc<1> GreensKernel::apply1d(Gaussian<1> &f) {
	double alpha, beta, ai, bi, ao, bo, p;
	const double *pos;

	ai = f.getExp();
	bi = f.getCoef();
	pos = f.getPos();

	alpha = kern.getExp(0);
	beta = kern.getCoef(0);

	p = ai + alpha;
	bo = bi * beta * sqrt(pi / p);
	ao = ai * alpha / p;

	return GaussFunc<1> (ao, bo, pos);
}

GaussExp<1> GreensKernel::apply1d(GaussExp<1> &f) {
	double alpha, beta, ai, bi, ao, bo, p;
	const double *pos;

	alpha = kern.getFunc(0).getExp();
	beta = kern.getFunc(0).getCoef();
	int nterms = f.size();
	GaussExp<1> convl(nterms);
	for (int i = 0; i < nterms; i++) {
		ai = f.getExp(i);
		bi = f.getCoef(i);
		pos = f.getPos(i);

		p = ai + alpha;
		bo = bi * beta * sqrt(pi / p);
		ao = ai * alpha / p;
		GaussFunc<1> g(ao, bo, pos);
		convl.setFunc(i, g);
	}
	return convl;
}

/** Convolution of a single Gaussian with a Kernel
 **/
GaussExp<1> GreensKernel::convolute(const Gaussian<1> &g) {
	const int kcomp = kern.size();

	double alpha, beta, ai, bi, ao, bo, p;
	const double *pos;

	GaussExp<1> convl(kcomp);
	ai = g.getExp();
	bi = g.getCoef();
	pos = g.getPos();

	for (int i = 0; i < kcomp; i++) {
		alpha = kern.getFunc(i).getExp();
		beta = kern.getFunc(i).getCoef();

		p = ai + alpha;

		bo = bi * beta * sqrt(pi / p);
		ao = ai * alpha / p;
		GaussFunc<1> f(ao, bo, pos);
		convl.setFunc(i, f);
	}

	return convl;
}

/** Convolution of a single Gaussian with the k:th Kernel component
 **/
GaussFunc<1> GreensKernel::convolute1(int k, const Gaussian<1> &g) {
	const int kcomp = kern.size();

	if (k < 0 or k > kcomp - 1) {
		return 0;
	}
	double alpha, beta, ai, bi, ao, bo, p;
	const double *pos;

	ai = g.getExp();
	bi = g.getCoef();
	pos = g.getPos();

	alpha = kern.getFunc(k).getExp();
	beta = kern.getFunc(k).getCoef();

	p = ai + alpha;
	bo = bi * beta * sqrt(pi / p);
	ao = ai * alpha / p;

	GaussFunc<1> f(ao, bo, pos);
	return f;
}

void GreensKernel::setCompMin(int n, bool reset) {
	if (n < 1 or n > compMax) {
		MSG_FATAL("compMin out of bounds!");
	}
	compMin = n;
	if (reset) {
		compKernRepr();
	}
}

void GreensKernel::setCompMax(int n, bool reset) {
	if (n < 1 or n < compMin) {
		MSG_FATAL("compMax out of bounds!");
	}
	compMax = n;
	if (reset) {
		compKernRepr();
	}
}

void GreensKernel::setRMin(double r, bool reset) {
	if (rMin > rMax) {
		MSG_FATAL("rMin > rMax");
	}
	rMin = r;
	if (reset) {
		compKernRepr();
	}
}

void GreensKernel::setRMax(double r, bool reset) {
	if (rMax < rMin) {
		MSG_FATAL("rMax < rMin");
	}
	rMax = r;
	if (reset) {
		compKernRepr();
	}
}

void GreensKernel::setExpPrec(double prec, bool reset) {
	expPrec = prec;
	if (reset) {
		compKernRepr();
	}
}

void GreensKernel::setProjPrec(double prec, bool reset) {
	projPrec = prec;
	if (reset) {
		compKernRepr();
	}
}

void GreensKernel::setKernelParam(double expPrec, double projPrec, double rMin,
double rMax, int nMin, int nMax) {
	if (expPrec > 0.0) {
		setExpPrec(expPrec, false);
	}
	if (projPrec > 0.0) {
		setProjPrec(projPrec, false);
	}
	if (rMin > 0.0) {
		setRMin(rMin, false);
	}
	if (rMax > 0.0) {
		setRMax(rMax, false);
	}
	if (nMin > 0) {
		setCompMin(nMin, false);
	}
	if (nMax > 0) {
		setCompMax(nMax, false);
	}
	compKernRepr();
}

void GreensKernel::setDefaultCompMin(int compMin) {
	if (compMin < 0) {
		MSG_ERROR("Requested compMin is negative: " << compMin);
		return;
	}
	defaultCompMin = compMin;
}

void GreensKernel::setDefaultCompMax(int compMax) {
	if (compMax < 0) {
		MSG_ERROR("Requested compMax is negative: " << compMax);
		return;
	}
	defaultCompMax = compMax;
}

void GreensKernel::setDefaultRMin(double rMin) {
	if (rMin < MIN_DIST) {
		MSG_ERROR("Requested rMin is too small: " << rMin);
		return;
	}
	defaultRMin = rMin;
}

void GreensKernel::setDefaultRMax(double rMax) {
	if (rMax > MAX_DIST) {
		MSG_ERROR("Requested rMax is too big: " << rMax);
		return;
	}
	defaultRMax = rMax;
}

void GreensKernel::setDefaultExpPrec(double expPrec) {
	if (expPrec < MIN_EPSILON) {
		MSG_ERROR("Requested expPrec is too small: " << expPrec);
		return;
	}
	if (expPrec > MAX_EPSILON) {
		MSG_ERROR("Requested expPrec is too big: " << expPrec);
		return;
	}
	defaultExpPrec = expPrec;
}

void GreensKernel::setDefaultProjPrec(double projPrec) {
	if (projPrec < MIN_EPSILON) {
		MSG_ERROR("Requested projPrec is too small: " << projPrec);
		return;
	}
	if (projPrec > MAX_EPSILON) {
		MSG_ERROR("Requested projPrec is too big: " << projPrec);
		return;
	}
	defaultProjPrec = projPrec;
}
