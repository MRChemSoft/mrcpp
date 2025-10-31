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

#include "GaussQuadrature.h"
#include "MRCPP/constants.h"
#include "MRCPP/macros.h"
#include "functions/LegendrePoly.h"
#include "utils/Printer.h"

using namespace Eigen;

namespace mrcpp {

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

void GaussQuadrature::rescaleRoots(VectorXd &rts, double a, double b, int inter) const {
    double transl = (b - a) / (double)inter;

    int k = 0;
    double pos = a;
    double xl = transl * 0.5;
    for (int i = 0; i < inter; i++) {
        for (int j = 0; j < this->order; j++) {
            rts(k) = this->unscaledRoots(j) * xl + pos + xl;
            ++k;
        }
        pos = pos + transl;
    }
}

void GaussQuadrature::rescaleWeights(VectorXd &wgts, double a, double b, int inter) const {
    double transl = (b - a) / (double)inter;

    int k = 0;
    double pos = a;
    double xl = transl * 0.5;
    for (int i = 0; i < inter; i++) {
        for (int j = 0; j < this->order; j++) {
            wgts(k) = this->unscaledWeights(j) * xl + pos + xl;
            ++k;
        }
        pos = pos + transl;
    }
}

void GaussQuadrature::calcScaledPtsWgts() {
    double transl = (this->B - this->A) / (double)this->intervals;

    int k = 0;
    double pos = this->A;
    double xl = transl * 0.5;
    for (int i = 0; i < this->intervals; i++) {
        for (int j = 0; j < this->order; j++) {
            this->roots(k) = this->unscaledRoots(j) * xl + pos + xl;
            this->weights(k) = this->unscaledWeights(j) * xl;
            ++k;
        }
        pos = pos + transl;
    }
}

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

    LegendrePoly legendrep(this->order, 1.0, 0.0);
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

double GaussQuadrature::integrate(RepresentableFunction<1> &func) const {
    double isum = 0.e0;
    Coord<1> r;
    for (int i = 0; i < this->npts; i++) {
        r[0] = this->roots(i);
        isum += this->weights(i) * func.evalf(r);
    }
    return isum;
}

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