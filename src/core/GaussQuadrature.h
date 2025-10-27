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

#pragma once

#include <Eigen/Core>

#include "functions/RepresentableFunction.h"

namespace mrcpp {

/**
 * @brief Maximum supported Gauss-Legendre order (per sub-interval).
 *
 * Implementation uses Newton iterations on Legendre polynomials in double
 * precision and is tuned for numerical stability up to this limit.
 */
const int MaxGaussOrder = 42;

/**
 * @brief Convergence tolerance for Newton's method when locating roots.
 */
static const double EPS = 3.0e-12;

/**
 * @brief Safety cap on Newton iterations per root.
 */
static const int NewtonMaxIter = 10;

/**
 * @brief Hard cap for a not-yet-implemented generic N-D integrator scaffold.
 * (Kept for legacy/planning; current code provides explicit 1D/2D/3D.)
 */
static const int MaxQuadratureDim = 7;

/**
 * @class GaussQuadrature
 * @brief Composite Gauss–Legendre quadrature on [A,B] with equal sub-intervals.
 *
 * What it represents
 * ------------------
 * A parameterized Gauss–Legendre rule over a (possibly partitioned) interval:
 *   - order      : number of Gauss nodes per sub-interval,
 *   - intervals  : number of equal pieces tiling [A,B],
 *   - roots      : all nodes over [A,B] for the composite rule (size npts),
 *   - weights    : corresponding weights (size npts).
 *
 * In addition, it stores the canonical (unscaled) Gauss nodes/weights on [-1,1]
 * so the rule can be remapped quickly if [A,B] or 'intervals' changes.
 *
 * Typical usage
 * -------------
 *   GaussQuadrature g(q, a, b, m);  // q = order, [a,b] bounds, m sub-intervals
 *   auto val1 = g.integrate(f1D);
 *   auto val2 = g.integrate(f2D);   // tensor-product rule (q*m in each axis)
 *
 * Notes
 * -----
 *  - “Composite” means we replicate the same order-q rule on each of the
 *    'intervals' equal sub-intervals, then sum the contributions.
 *  - setBounds() / setIntervals() preserve the canonical [-1,1] rule and
 *    rebuild the scaled (roots,weights) for the new configuration.
 */
class GaussQuadrature final {
public:
    /**
     * @brief Construct a Gauss–Legendre quadrature rule.
     * @param k      Order (nodes per sub-interval), 0 ≤ k ≤ MaxGaussOrder.
     * @param a      Lower bound A (default -1).
     * @param b      Upper bound B (default  1).
     * @param inter  Number of equal sub-intervals (default 1, must be ≥ 1).
     *
     * Effects (see .cpp):
     *  - Builds canonical nodes/weights on [-1,1] via Newton’s method.
     *  - Scales/replicates them to fill (roots,weights) over [A,B].
     */
    GaussQuadrature(int k, double a = -1.0, double b = 1.0, int inter = 1);

    /**
     * @name Tensor-product integration helpers
     * @{
     * @brief Integrate a RepresentableFunction using the prepared rule.
     *
     * 1D:   ∑_i w_i f(x_i)
     * 2D:   ∑_i ∑_j w_i w_j f(x_i, x_j)
     * 3D:   ∑_i ∑_j ∑_k w_i w_j w_k f(x_i, x_j, x_k)
     */
    double integrate(RepresentableFunction<1> &func) const;
    double integrate(RepresentableFunction<2> &func) const;
    double integrate(RepresentableFunction<3> &func) const;
    /** @} */

    /**
     * @brief Set the number of equal sub-intervals and rebuild (roots,weights).
     * @param i New number of sub-intervals (≥ 1).
     *
     * Reallocates global arrays to size npts = order * intervals and remaps
     * the canonical [-1,1] rule accordingly.
     */
    void setIntervals(int i);

    /**
     * @brief Set integration bounds [a,b] and rebuild (roots,weights).
     * @param a Lower bound
     * @param b Upper bound (must satisfy a < b)
     */
    void setBounds(double a, double b);

    /** @return Number of sub-intervals tiling [A,B]. */
    int getIntervals() const { return this->intervals; }
    /** @return Upper bound B. */
    double getUpperBound() const { return this->B; }
    /** @return Lower bound A. */
    double getLowerBound() const { return this->A; }

    /** @return Composite-rule nodes over [A,B] (size npts). */
    const Eigen::VectorXd &getRoots() const { return this->roots; }
    /** @return Composite-rule weights over [A,B] (size npts). */
    const Eigen::VectorXd &getWeights() const { return this->weights; }

    /** @return Canonical Gauss nodes on [-1,1] (size order). */
    const Eigen::VectorXd &getUnscaledRoots() const { return this->unscaledRoots; }
    /** @return Canonical Gauss weights on [-1,1] (size order). */
    const Eigen::VectorXd &getUnscaledWeights() const { return this->unscaledWeights; }

protected:
    // ---- Parameters describing the rule ----
    int order;        ///< Nodes per sub-interval (q)
    double A;         ///< Lower integration bound
    double B;         ///< Upper integration bound
    int intervals;    ///< Number of equal sub-intervals tiling [A,B]
    int npts;         ///< Total nodes = order * intervals

    // ---- Scaled (composite) rule on [A,B] ----
    Eigen::VectorXd roots;          ///< All nodes over [A,B] (size npts)
    Eigen::VectorXd weights;        ///< All weights over [A,B] (size npts)

    // ---- Canonical rule on [-1,1] ----
    Eigen::VectorXd unscaledRoots;   ///< Nodes on [-1,1] (size 'order')
    Eigen::VectorXd unscaledWeights; ///< Weights on [-1,1] (size 'order')

    /**
     * @brief Map canonical nodes onto [a,b] replicated over @p inter blocks.
     * @param rts   Output vector of length inter*order.
     * @param a,b   Interval bounds.
     * @param inter Number of sub-intervals (default 1).
     *
     * Each block is an affine image of [-1,1] with width (b-a)/inter.
     */
    void rescaleRoots(Eigen::VectorXd &rts, double a, double b, int inter = 1) const;

    /**
     * @brief Map canonical weights onto [a,b] replicated over @p inter blocks.
     * @param wgts  Output vector of length inter*order.
     * @param a,b   Interval bounds.
     * @param inter Number of sub-intervals (default 1).
     *
     * Weights scale by the Jacobian of the affine mapping (factor 0.5*(b-a)/inter).
     */
    void rescaleWeights(Eigen::VectorXd &wgts, double a, double b, int inter = 1) const;

    /**
     * @brief Rebuild (roots,weights) on [A,B] for the current 'intervals'.
     *
     * Uses the stored canonical (unscaled) rule and performs replication.
     */
    void calcScaledPtsWgts();

    /**
     * @brief Compute canonical Gauss–Legendre nodes/weights on [-1,1].
     * @return 1 on success; 0 if Newton iteration failed to converge.
     *
     * Uses Newton’s method on Legendre polynomials with symmetric pairing:
     * computes half the roots in (0,1) and reflects them about 0.
     */
    int calcGaussPtsWgts();

    /**
     * @brief Planned recursive N-D integration (not implemented).
     * @return No return; aborts at runtime if called.
     */
    double integrate_nd(RepresentableFunction<3> &func, int axis = 0) const;
};

} // namespace mrcpp