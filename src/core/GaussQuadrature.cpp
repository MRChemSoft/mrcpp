/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
 */

#include "GaussQuadrature.h"
#include "MRCPP/constants.h"
#include "functions/LegendrePoly.h"
#include "MRCPP/macros.h"
#include "utils/Printer.h"

using namespace Eigen;

namespace mrcpp {

/** Constructor for Gauss-Legendre quadrature.
 *
 * \param order Polynominal order
 * \param a Lower bound of validity
 * \param b Upper bound of validity
 * \param intervals Number of intervals to divde |a-b| into
 */
GaussQuadrature::GaussQuadrature(int k, double a, double b, int inter) {
    this->order = k;
    this->A = a;
    this->B = b;
    this->intervals = inter;

    if (this->order < 0 || this->order > MaxGaussOrder) {
        MSG_ERROR("Gauss quadrature order " << this->order << " is larger than the maximum of " << MaxGaussOrder);
    }
    if (a >= b) { MSG_ERROR("Invalid Gauss interval, a > b."); }
    if (this->intervals < 1) { MSG_ERROR("Invalid number of intervals, intervals < 1"); }
    this->npts = this->order * this->intervals;
    this->roots = VectorXd::Zero(this->npts);
    this->weights = VectorXd::Zero(this->npts);
    this->unscaledRoots = VectorXd::Zero(this->order);
    this->unscaledWeights = VectorXd::Zero(this->order);
    // set up unscaled Gauss points and weights ( interval ]-1,1[)
    if (calcGaussPtsWgts() != 1) { MSG_ERROR("Setup of Gauss-Legendre weights failed.") }
    calcScaledPtsWgts();
}

void GaussQuadrature::setBounds(double a, double b) {
    if (std::abs(this->A - a) < MachineZero and std::abs(this->B - b) < MachineZero) { return; }
    if (a >= b) { MSG_ERROR("Invalid bounds: a > b"); }
    this->A = a;
    this->B = b;
    calcScaledPtsWgts();
}

void GaussQuadrature::setIntervals(int i) {
    if (i == this->intervals) { return; }
    if (i < 1) { MSG_ERROR("Invalid number of integration intervals: " << i); }
    this->intervals = i;
    this->npts = this->order * this->intervals;
    this->roots = VectorXd::Zero(this->npts);
    this->weights = VectorXd::Zero(this->npts);
    calcScaledPtsWgts();
}

/** Calculate scaled distribution of roots for Gauss-Legendre
 * quadrature on on ]a,b[. The number of quadrature points on the interval
 * is scale*(order+1).
 */
void GaussQuadrature::rescaleRoots(VectorXd &rts, double a, double b, int inter) const {
    // length of one block
    double transl = (b - a) / (double)inter;

    int k = 0;
    double pos = a;
    double xl = transl * 0.5;
    // scale and translate Gauss points and weights
    for (int i = 0; i < inter; i++) {
        for (int j = 0; j < this->order; j++) {
            rts(k) = this->unscaledRoots(j) * xl + pos + xl;
            ++k;
        }
        pos = pos + transl;
    }
}

/** Calculate scaled distribution of weights for Gauss-Legendre
 * quadrature on on ]a,b[. The number of quadrature points on the interval
 * is scale*(order+1).
 */
void GaussQuadrature::rescaleWeights(VectorXd &wgts, double a, double b, int inter) const {
    // length of one block
    double transl = (b - a) / (double)inter;

    int k = 0;
    double pos = a;
    double xl = transl * 0.5;
    // scale and translate Gauss points and weights
    for (int i = 0; i < inter; i++) {
        for (int j = 0; j < this->order; j++) {
            wgts(k) = this->unscaledWeights(j) * xl + pos + xl;
            ++k;
        }
        pos = pos + transl;
    }
}

/** Calculate scaled distribution of points and weights for Gauss-Legendre
 * quadrature on on ]a,b[. The number of quadrature points on the interval
 * is scale*(order+1).
 */
void GaussQuadrature::calcScaledPtsWgts() {
    // length of one block
    double transl = (this->B - this->A) / (double)this->intervals;

    int k = 0;
    double pos = this->A;
    double xl = transl * 0.5;
    // scale and translate Gauss points and weights
    for (int i = 0; i < this->intervals; i++) {
        for (int j = 0; j < this->order; j++) {
            this->roots(k) = this->unscaledRoots(j) * xl + pos + xl;
            this->weights(k) = this->unscaledWeights(j) * xl;
            ++k;
        }
        pos = pos + transl;
    }
}

/** Calulate distribution of points and weights for Guass-Legendre quadrature on
 * ]-1,1[.
 *
 * Find quadrature points and weights by solving for the roots of
 * Legendre polynomials using Newtons method. Using double precison the
 * maximum stable order is currently set to 13. Return 1 on success, 0 on failure.
 *
 */
int GaussQuadrature::calcGaussPtsWgts() {
    int K;
    if (this->order % 2 == 0) {
        K = this->order / 2;
    } else {
        K = (this->order + 1) / 2;
    }

    double a = -1.0;
    double b = 1.0;

    double xm = (b + a) * 0.5;
    double xl = (b - a) * 0.5;

    LegendrePoly legendrep(this->order, 1.0, 0.0); // Interval [-1,1]
    Vector2d lp;

    for (int i = 0; i < K; i++) {
        double z = cos(pi * (i + 0.75) / (this->order + 0.5));
        int iter;
        for (iter = 0; iter < NewtonMaxIter; iter++) {
            lp = legendrep.firstDerivative(z);

            double z1 = z;
            z = z1 - lp(0) / lp(1);
            if (std::abs(z - z1) <= EPS) { break; }
        }
        if (iter == NewtonMaxIter) { return 0; }

        this->unscaledRoots(i) = xm - xl * z;
        this->unscaledRoots(order - 1 - i) = xm + xl * z;

        this->unscaledWeights(i) = 2.e0 * xl / ((1.e0 - z * z) * lp(1) * lp(1));
        this->unscaledWeights(order - 1 - i) = this->unscaledWeights(i);
    }
    return 1;
}

/** Integrate a 1D-function f(x) using quadrature */
double GaussQuadrature::integrate(RepresentableFunction<1> &func) const {
    double isum = 0.e0;
    Coord<1> r;
    for (int i = 0; i < this->npts; i++) {
        r[0] = this->roots(i);
        isum += this->weights(i) * func.evalf(r);
    }
    return isum;
}

/** Integrate a 2D-function f(x1, x2) using quadrature */
double GaussQuadrature::integrate(RepresentableFunction<2> &func) const {
    Coord<2> r;
    double isum = 0.e0;
    for (int i = 0; i < this->npts; i++) {
        double jsum = 0.e0;
        r[0] = this->roots(i);
        for (int j = 0; j < this->npts; j++) {
            r[1] = this->roots(j);
            jsum += this->weights(j) * func.evalf(r);
        }
        isum += jsum * this->weights(i);
    }
    return isum;
}

/** Integrate a 3D-function f(x1, x2, x3) using quadrature */
double GaussQuadrature::integrate(RepresentableFunction<3> &func) const {
    Coord<3> r;

    double isum = 0.e0;
    for (int i = 0; i < this->npts; i++) {
        double jsum = 0.e0;
        r[0] = this->roots(i);
        for (int j = 0; j < this->npts; j++) {
            double ksum = 0.e0;
            r[1] = this->roots(j);
            for (int k = 0; k < this->npts; k++) {
                r[2] = this->roots(k);
                ksum += this->weights(k) * func.evalf(r);
            }
            jsum += ksum * this->weights(j);
        }
        isum += jsum * this->weights(i);
    }
    return isum;
}

/** Integrate a ND-function f(x1,...), allowing for different
 * quadrature in each dimension.
 *
 * This function has been implemented using a recursive algorithm.
 */
double GaussQuadrature::integrate_nd(RepresentableFunction<3> &func, int axis) const {
    NOT_IMPLEMENTED_ABORT;
    NEEDS_TESTING

    /*
    double sum;
    static double r[MaxQuadratureDim];

    sum = 0.e0;
    for (int i = 0; i < this->npts; i++) {
        r[axis] = this->roots(i);
        if (axis < 2) {
            sum += integrate_nd(func, axis) * this->weights(i);
        } else {
            sum += this->weights(i) * func.evalf(r);
        }
    }
    return sum;
    */
}

} // namespace mrcpp
