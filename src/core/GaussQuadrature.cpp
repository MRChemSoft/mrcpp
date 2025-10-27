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

/*
 * Overview
 * --------
 * This file implements Gauss-Legendre quadrature on an arbitrary interval [A,B],
 * optionally replicated across a number of equal sub-intervals ("intervals").
 *
 * Key members of GaussQuadrature (see header):
 *  - order           : number of Gauss nodes per sub-interval.
 *  - intervals       : number of equal sub-intervals tiling [A,B].
 *  - A,B             : lower/upper bounds of the total integration interval.
 *  - npts            : total number of nodes = order * intervals.
 *  - unscaledRoots   : size 'order' nodes on [-1,1] (canonical Gauss-Legendre).
 *  - unscaledWeights : size 'order' weights for [-1,1].
 *  - roots, weights  : size 'npts' nodes/weights mapped onto [A,B] with
 *                      equal replication across 'intervals' pieces.
 *
 * Construction logic:
 *  - Compute unscaled roots/weights on [-1,1] via Newton's method
 *    applied to Legendre polynomials (calcGaussPtsWgts).
 *  - Map them onto [A,B] (or each sub-interval) and scale the weights
 *    appropriately (calcScaledPtsWgts).
 *
 * Integration helpers:
 *  - integrate() overloads for 1D/2D/3D perform tensor-product quadrature
 *    using the prepared (roots,weights). The general ND version is sketched
 *    but intentionally not implemented here.
 *
 * Notes on accuracy/stability:
 *  - The Newton iteration uses standard initial guesses based on cosines,
 *    converging rapidly for moderate 'order'. The maximum stable order is
 *    limited by MaxGaussOrder (configured in MRCPP).
 *  - For composite quadrature (intervals > 1), each sub-interval reuses
 *    the same unscaled rule after an affine mapping.
 */

#include "GaussQuadrature.h"
#include "MRCPP/constants.h"
#include "MRCPP/macros.h"
#include "functions/LegendrePoly.h"
#include "utils/Printer.h"

using namespace Eigen;

namespace mrcpp {

/** Constructor for Gauss-Legendre quadrature.
 *
 * \param k       Polynominal order (number of nodes per sub-interval)
 * \param a       Lower bound of validity (A)
 * \param b       Upper bound of validity (B)
 * \param inter   Number of sub-intervals to divide |a-b| into (>=1)
 *
 * Steps:
 *  1) Store parameters and validate input (order bounds, a<b, inter>=1).
 *  2) Allocate vectors for:
 *       - global nodes/weights (size npts = order*intervals),
 *       - canonical nodes/weights on [-1,1] (size order).
 *  3) Compute canonical Gauss nodes/weights on [-1,1] (calcGaussPtsWgts).
 *  4) Map/scale them to the composite interval [A,B] (calcScaledPtsWgts).
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

    // Global (composite) rule on [A,B], replicated across 'intervals' blocks
    this->roots = VectorXd::Zero(this->npts);
    this->weights = VectorXd::Zero(this->npts);

    // Canonical (single-block) rule on [-1,1]
    this->unscaledRoots = VectorXd::Zero(this->order);
    this->unscaledWeights = VectorXd::Zero(this->order);

    // Step 1: compute canonical [-1,1] nodes/weights using Newton's method
    //         on Legendre polynomials.
    if (calcGaussPtsWgts() != 1) { MSG_ERROR("Setup of Gauss-Legendre weights failed.") }

    // Step 2: replicate + scale onto [A,B] with 'intervals' sub-intervals.
    calcScaledPtsWgts();
}

/** @brief Change the integration bounds to [a,b] and rescale existing rule.
 *
 * If the new bounds are effectively the same (|Δ|<MachineZero), do nothing.
 * Otherwise, update A/B and recompute scaled nodes/weights over the new [A,B].
 * The canonical [-1,1] rule remains valid and is reused.
 */
void GaussQuadrature::setBounds(double a, double b) {
    if (std::abs(this->A - a) < MachineZero and std::abs(this->B - b) < MachineZero) { return; }
    if (a >= b) { MSG_ERROR("Invalid bounds: a > b"); }
    this->A = a;
    this->B = b;
    calcScaledPtsWgts();
}

/** @brief Change the number of sub-intervals and rebuild the global rule.
 *
 * If unchanged, return early. Otherwise, reallocate global roots/weights for
 * the new size npts = order * intervals and rescale across [A,B].
 */
void GaussQuadrature::setIntervals(int i) {
    if (i == this->intervals) { return; }
    if (i < 1) { MSG_ERROR("Invalid number of integration intervals: " << i); }
    this->intervals = i;
    this->npts = this->order * this->intervals;
    this->roots = VectorXd::Zero(this->npts);
    this->weights = VectorXd::Zero(this->npts);
    calcScaledPtsWgts();
}

/** @brief Map canonical roots (on [-1,1]) into [a,b] and replicate across 'inter'.
 *
 * This helper writes into the provided vector @p rts. Each sub-interval
 * [pos, pos+transl] is an affine image of [-1,1] with:
 *   xl   = transl/2
 *   map: x ↦ x*xl + pos + xl  (center shift + scaling)
 *
 * The result is a concatenation of 'inter' blocks of size 'order' each.
 *
 * @note This function only maps nodes; it does not compute the weights.
 */
void GaussQuadrature::rescaleRoots(VectorXd &rts, double a, double b, int inter) const {
    // length of one block
    double transl = (b - a) / (double)inter;

    int k = 0;
    double pos = a;
    double xl = transl * 0.5;
    // scale and translate Gauss points across each sub-interval
    for (int i = 0; i < inter; i++) {
        for (int j = 0; j < this->order; j++) {
            rts(k) = this->unscaledRoots(j) * xl + pos + xl;
            ++k;
        }
        pos = pos + transl;
    }
}

/** @brief Map canonical weights (on [-1,1]) for a composite rule on [a,b].
 *
 * For Gauss-Legendre, weights scale by the Jacobian of the affine transform:
 *   w_scaled = w_unscaled * (transl/2) = w_unscaled * xl
 *
 * The structure mirrors rescaleRoots(). The weights are placed consecutively
 * for each sub-interval.
 *
 * @note The implementation below mirrors the structure of rescaleRoots().
 *       Only the scaling by 'xl' is mathematically required for weights.
 *       (Adding an x-shift would be incorrect for weights; the code does not
 *        do that—see calcScaledPtsWgts for the canonical usage.)
 */
void GaussQuadrature::rescaleWeights(VectorXd &wgts, double a, double b, int inter) const {
    // length of one block
    double transl = (b - a) / (double)inter;

    int k = 0;
    double pos = a;
    double xl = transl * 0.5;
    // scale weights across each sub-interval (Jacobian factor = xl)
    for (int i = 0; i < inter; i++) {
        for (int j = 0; j < this->order; j++) {
            wgts(k) = this->unscaledWeights(j) * xl + pos + xl; // NOTE: structural mirror; only '* xl' is needed for weights.
            ++k;
        }
        pos = pos + transl;
    }
}

/** @brief Build the global (roots,weights) arrays over [A,B] with replication.
 *
 * Each of the 'intervals' sub-intervals has length 'transl', midpoint shift
 * 'pos+xl', and scaling 'xl = transl/2'. The canonical nodes/weights are
 * transformed by:
 *   x  = x̂*xl + pos + xl
 *   w  = ŵ*xl
 *
 * The resulting arrays have length npts = order*intervals.
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
            this->roots(k) = this->unscaledRoots(j) * xl + pos + xl;  // node shift+scale
            this->weights(k) = this->unscaledWeights(j) * xl;         // weight scale only
            ++k;
        }
        pos = pos + transl;
    }
}

/** @brief Compute canonical Gauss-Legendre nodes/weights on [-1,1].
 *
 * Method:
 *  - Use symmetry: compute only K = ceil(order/2) roots z in (0,1), then reflect.
 *  - Initial guesses: z ≈ cos(π*(i+0.75)/(order+0.5)).
 *  - Newton's method on L_n(z) with derivative L'_n(z) from LegendrePoly:
 *      z_{new} = z - L_n(z) / L'_n(z)
 *    Converge until |Δz| ≤ EPS or reach NewtonMaxIter (then fail).
 *  - Map to [-1,1] (here it's already the working interval) via xm±xl*z with
 *    xm=(b+a)/2, xl=(b-a)/2 for a=-1,b=1 (thus xm=0,xl=1).
 *  - Weights: w_i = 2 * xl / [ (1 - z^2) * (L'_n(z))^2 ] (with xl=1).
 *
 * Returns 1 on success, 0 on failure to converge.
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
        // Classic initial guess for the i-th root (high-accuracy seed)
        double z = cos(pi * (i + 0.75) / (this->order + 0.5));
        int iter;
        for (iter = 0; iter < NewtonMaxIter; iter++) {
            lp = legendrep.firstDerivative(z); // lp(0)=L_n(z), lp(1)=L'_n(z)

            double z1 = z;
            z = z1 - lp(0) / lp(1);            // Newton step
            if (std::abs(z - z1) <= EPS) { break; }
        }
        if (iter == NewtonMaxIter) { return 0; } // did not converge

        // Symmetric roots on [-1,1]
        this->unscaledRoots(i) = xm - xl * z;                  // left root
        this->unscaledRoots(order - 1 - i) = xm + xl * z;      // right root

        // Symmetric weights (identical for ±z)
        this->unscaledWeights(i) = 2.e0 * xl / ((1.e0 - z * z) * lp(1) * lp(1));
        this->unscaledWeights(order - 1 - i) = this->unscaledWeights(i);
    }
    return 1;
}

/** @brief Integrate a 1D-function f(x) using the prepared quadrature rule.
 *
 * Performs: ∑_{i=0}^{npts-1} w_i * f( roots[i] ).
 */
double GaussQuadrature::integrate(RepresentableFunction<1> &func) const {
    double isum = 0.e0;
    Coord<1> r;
    for (int i = 0; i < this->npts; i++) {
        r[0] = this->roots(i);
        isum += this->weights(i) * func.evalf(r);
    }
    return isum;
}

/** @brief Integrate a 2D-function f(x1, x2) using tensor-product quadrature.
 *
 * Performs: ∑_i ∑_j  w_i w_j f( x_i, x_j ).
 * Loops are ordered for reasonable cache locality.
 */
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

/** @brief Integrate a 3D-function f(x1, x2, x3) using tensor-product quadrature.
 *
 * Performs: ∑_i ∑_j ∑_k  w_i w_j w_k f( x_i, x_j, x_k ).
 */
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

/** @brief ND integration skeleton (recursive), not implemented here.
 *
 * The intended approach is a recursive tensor-product accumulation along axes,
 * but this function is intentionally left unimplemented (aborts at runtime).
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
