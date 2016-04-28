/*
 * 
 *
 *  \date Jul 7, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include <cmath>
#include "HelmholtzKernel.h"
#include "macros.h"

using namespace std;

/** generate an approximation of the 3d helmholtz kernel expanded in gaussian functions
 */
void HelmholtzKernel::compKernRepr() {
    double r0 = this->rMin;
    double r1 = this->rMax;
    double mu_tilde = mu*r1;
    double epsilon = this->expPrec;

    // Set the truncation limits s1,s2 of the integral (integrate over [s1,s2])
    // for achieving relative error epsilon
    double t = max((-2.5L * log(epsilon)), 5.0L);
    double s1 = -log(4 * t / (mu_tilde * mu_tilde)) / 2;
    double s2 = log(t / (r0 * r0)) / 2;

    int p = 0;
    if (delta) {
	NOT_IMPLEMENTED_ABORT
	p = 1;
    }

    // Now, set the proper step size h for use in the trapezoidal rule
    // for given MU
    double h = 1.0 / (0.20L - 0.47L * log10(epsilon));
    int n_exp = (int) ceil((s2 - s1) / h) + 1;
    if (n_exp > MAX_SEP_RANK)
	MSG_FATAL("Maximum separation rank exceeded.");

    kern = GaussExp<1> (n_exp);
    for (int i = 0; i < n_exp; i++) {
    	double arg = s1 + h * i;
    	double temp = -arg * 2.0;
	//cout << "i " << i << " arg " << arg << " temp " << temp << endl;
	double temp2 = -mu_tilde * mu_tilde * exp(temp) / 4.0 + arg;
	double beta = deltaMu * (h * (2.0 / root_pi) * exp(temp2));
	double temp3 = 2.0L * arg;
	double alpha = exp(temp3);
	double pos[1] = {0.0};
	int power[1] = {p};

	alpha *= 1.0/(r1*r1);
	beta *= 1.0/r1;

	GaussFunc<1> gFunc(alpha, beta, pos, power);
	kern.setFunc(i, gFunc);
    }
    double beta = kern.getCoef(0)/2.0;
    kern.setCoef(0, beta);
    beta = kern.getCoef(n_exp - 1)/2.0;
    kern.setCoef(n_exp - 1, beta);
    kern.calcSquareNorm();
}
