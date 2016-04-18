/*
 * 
 *
 *  \date Jul 7, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "PoissonKernel.h"

#include <cmath>
#include <cassert>

/** generate an approximation of the 3d poisson kernel expanded in
 * gaussian functions this routine assumes that the expansion be centered
 */
void PoissonKernel::compKernRepr() {
    double r0 = this->rMin;
    double r1 = this->rMax;
    double epsilon = this->expPrec;

    double t1 = 1.0L;
    while ((2 * t1 * exp(-t1)) > epsilon) {
	t1 *= 1.1L;
    }
    double t2 = 1.0L;
    while ((sqrt(t2) * exp(-t2) / r0) > epsilon) {
	t2 *= 1.1L;
    }

    // Set the truncation limits s1,s2 of the integral (integrate over [s1,s2])
    // for achieving relative error epsilon
    double s1 = -log(2 * t1);
    double s2 = log(t2 / (r0 * r0)) / 2;

    // Now, set the step size h for use in the trapezoidal rule for given MU
    double h = 1 / (0.2L - 0.47L * log10(epsilon));
    int n_exp = (int) ceil((s2 - s1) / h) + 1;

    if (n_exp > MAX_SEP_RANK) {
    	MSG_FATAL("Maximum separation rank exceeded.");
    }

    kern = GaussExp<1> (n_exp);
    for (int i = 0; i < n_exp; i++) {
	double arg = s1 + h * i;
	double sinharg = sinh(arg);
	double cosharg = cosh(arg);
	double onepexp = 1.0 + exp(-sinharg);

    	double alpha = 4.0L * (sinharg+log(onepexp)) * (sinharg+log(onepexp));
	double beta = h * (4.0L / root_pi) * cosharg / onepexp;
	double pos[1] = {0.0};
	int power[1] = {0};

	alpha *= 1.0/(r1*r1);
	beta *= 1.0/r1;

	GaussFunc<1> g(alpha, beta, pos, power);
	kern.setFunc(i, g);
    }

    double beta = kern.getCoef(0) / 2.0;
    kern.setCoef(0, beta);
    beta = kern.getCoef(n_exp - 1) / 2.0;
    kern.setCoef(n_exp - 1, beta);
    kern.calcSquareNorm();
}


