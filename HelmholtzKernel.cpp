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
#include "GaussFunc.h"

using namespace std;

/** generate an approximation of the 3d helmholtz kernel expanded in gaussian functions
 */
void HelmholtzKernel::initializeKernel() {
    double r0 = this->rMin;
    double r1 = this->rMax;
    double mu_tilde = this->mu*r1;

    // Set the truncation limits s1,s2 of the integral (integrate over [s1,s2])
    // for achieving relative error epsilon
    double t = max((-2.5L * log(this->epsilon)), 5.0L);
    double s1 = -log(4 * t / (mu_tilde * mu_tilde)) / 2;
    double s2 = log(t / (r0 * r0)) / 2;

    // Now, set the proper step size h for use in the trapezoidal rule
    // for given MU
    double h = 1.0 / (0.20L - 0.47L * log10(this->epsilon));
    int n_exp = (int) ceil((s2 - s1) / h) + 1;
    if (n_exp > MaxSepRank) MSG_FATAL("Maximum separation rank exceeded.");

    for (int i = 0; i < n_exp; i++) {
        double arg = s1 + h * i;
        double temp = -arg * 2.0;
        double temp2 = -mu_tilde * mu_tilde * exp(temp) / 4.0 + arg;
        double beta = (h * (2.0 / root_pi) * exp(temp2));
        double temp3 = 2.0L * arg;
        double alpha = exp(temp3);
        double pos = 0.0;

        alpha *= 1.0/(r1*r1);
        if (i == 0 or i == (n_exp - 1)) {
            beta *= 1.0/(2.0*r1);
        } else {
            beta *= 1.0/r1;
        }
        beta = pow(beta, 1.0/3.0);

        GaussFunc<1> gFunc(alpha, beta, &pos);
        this->append(gFunc);
    }
    this->calcSquareNorm();
}
