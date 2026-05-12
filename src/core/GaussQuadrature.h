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

/** @brief Maximum supported Gauss-Legendre order per sub-interval */
const int MaxGaussOrder = 42;

/** @brief Convergence tolerance for Newton's method when locating Gauss roots */
static const double EPS = 3.0e-12;

/** @brief Maximum number of Newton iterations allowed per root */
static const int NewtonMaxIter = 10;

/** @brief Dimension cap reserved for a planned N-D integrator (not yet implemented) */
static const int MaxQuadratureDim = 7;

/**
 * @class GaussQuadrature
 * @brief Composite Gauss–Legendre quadrature rule on \f$ [A,B] \f$ divided into equal sub-intervals
 *
 * @details
 * The rule stores @p order Gauss nodes and weights per sub-interval and replicates them over @p intervals
 * equal pieces of \f$ [A,B] \f$, giving a total of \f$ \text{npts} = \text{order} \times \text{intervals} \f$
 * nodes. The canonical (unscaled) nodes and weights on \f$ [-1,1] \f$ are kept internally so that the
 * composite rule can be rebuilt quickly when the bounds or interval count change via setBounds() or
 * setIntervals(). Integration of RepresentableFunction objects in 1D, 2D, and 3D is provided as
 * tensor-product sums.
 *
 * @see QuadratureCache for the singleton that caches GaussQuadrature objects by order
 */
class GaussQuadrature final {
public:
    /**
     * @brief Construct a composite Gauss–Legendre quadrature rule
     * @param k     Number of nodes per sub-interval (\f$ 0 \leq k \leq \f$ MaxGaussOrder)
     * @param a     Lower integration bound (default \f$ -1 \f$)
     * @param b     Upper integration bound (default \f$ 1 \f$; must satisfy \f$ a < b \f$)
     * @param inter Number of equal sub-intervals (default 1; must be \f$ \geq 1 \f$)
     *
     * @details Computes canonical nodes and weights on \f$ [-1,1] \f$ via Newton’s method on Legendre
     * polynomials (using symmetric pairing to halve the work), then scales and replicates them over
     * the @p inter sub-intervals of \f$ [a,b] \f$
     */
    GaussQuadrature(int k, double a = -1.0, double b = 1.0, int inter = 1);

    /**
     * @brief Integrate a 1D RepresentableFunction: \f$ \sum_i w_i f(x_i) \f$
     * @param[in] func The function to integrate
     * @return Numerical integral
     */
    double integrate(RepresentableFunction<1> &func) const;
    /**
     * @brief Integrate a 2D RepresentableFunction via tensor-product rule:
     *        \f$ \sum_i \sum_j w_i w_j f(x_i, x_j) \f$
     * @param[in] func The function to integrate
     * @return Numerical integral
     */
    double integrate(RepresentableFunction<2> &func) const;
    /**
     * @brief Integrate a 3D RepresentableFunction via tensor-product rule:
     *        \f$ \sum_i \sum_j \sum_k w_i w_j w_k f(x_i, x_j, x_k) \f$
     * @param[in] func The function to integrate
     * @return Numerical integral
     */
    double integrate(RepresentableFunction<3> &func) const;

    /**
     * @brief Update the number of sub-intervals and rebuild the composite rule
     * @param i New number of sub-intervals (\f$ \geq 1 \f$)
     *
     * @details Resizes #roots and #weights to \f$ \text{order} \times i \f$ and recomputes
     * the scaled nodes and weights over \f$ [A,B] \f$
     */
    void setIntervals(int i);

    /**
     * @brief Update the integration bounds and rebuild the composite rule
     * @param a New lower bound (must satisfy \f$ a < b \f$)
     * @param b New upper bound
     */
    void setBounds(double a, double b);

    /** @return Number of sub-intervals tiling [A,B] */
    int getIntervals() const { return this->intervals; }
    /** @return Upper bound B */
    double getUpperBound() const { return this->B; }
    /** @return Lower bound A */
    double getLowerBound() const { return this->A; }

    /** @return Composite-rule nodes over [A,B] (size npts) */
    const Eigen::VectorXd &getRoots() const { return this->roots; }
    /** @return Composite-rule weights over [A,B] (size npts) */
    const Eigen::VectorXd &getWeights() const { return this->weights; }

    /** @return Canonical Gauss nodes on [-1,1] (size order) */
    const Eigen::VectorXd &getUnscaledRoots() const { return this->unscaledRoots; }
    /** @return Canonical Gauss weights on [-1,1] (size order) */
    const Eigen::VectorXd &getUnscaledWeights() const { return this->unscaledWeights; }

protected:
    int order;     ///< Number of Gauss nodes per sub-interval
    double A;      ///< Lower integration bound
    double B;      ///< Upper integration bound
    int intervals; ///< Number of equal sub-intervals tiling \f$ [A,B] \f$
    int npts;      ///< Total number of nodes: \f$ \text{order} \times \text{intervals} \f$

    Eigen::VectorXd roots;          ///< All composite-rule nodes over \f$ [A,B] \f$ (size #npts)
    Eigen::VectorXd weights;        ///< All composite-rule weights over \f$ [A,B] \f$ (size #npts)

    Eigen::VectorXd unscaledRoots;   ///< Canonical Gauss nodes on \f$ [-1,1] \f$ (size #order)
    Eigen::VectorXd unscaledWeights; ///< Canonical Gauss weights on \f$ [-1,1] \f$ (size #order)

    /**
     * @brief Map the canonical nodes from \f$ [-1,1] \f$ onto @p inter equal blocks of \f$ [a,b] \f$
     * @param[out] rts  Output vector of length \f$ \text{inter} \times \text{order} \f$
     * @param      a    Lower bound of the target interval
     * @param      b    Upper bound of the target interval
     * @param      inter Number of sub-intervals (default 1)
     */
    void rescaleRoots(Eigen::VectorXd &rts, double a, double b, int inter = 1) const;

    /**
     * @brief Map the canonical weights from \f$ [-1,1] \f$ onto @p inter equal blocks of \f$ [a,b] \f$
     * @param[out] wgts Output vector of length \f$ \text{inter} \times \text{order} \f$
     * @param      a    Lower bound of the target interval
     * @param      b    Upper bound of the target interval
     * @param      inter Number of sub-intervals (default 1)
     *
     * @details Each weight is scaled by the Jacobian \f$ \tfrac{b-a}{2\,\text{inter}} \f$
     */
    void rescaleWeights(Eigen::VectorXd &wgts, double a, double b, int inter = 1) const;

    /**
     * @brief Rebuild #roots and #weights on \f$ [A,B] \f$ using the current #intervals value
     *
     * @details Applies rescaleRoots() and rescaleWeights() using the stored canonical rule
     */
    void calcScaledPtsWgts();

    /**
     * @brief Compute canonical Gauss–Legendre nodes and weights on \f$ [-1,1] \f$
     * @return 1 on success; 0 if Newton iteration failed to converge for any root
     *
     * @details Uses Newton’s method on Legendre polynomials. Exploits symmetry about \f$ x=0 \f$
     * to compute only half the nodes, then reflects them.
     */
    int calcGaussPtsWgts();

    /**
     * @brief Placeholder for a future recursive N-D integrator
     * @return This method aborts at runtime; it is not yet implemented
     */
    double integrate_nd(RepresentableFunction<3> &func, int axis = 0) const;
};

} // namespace mrcpp