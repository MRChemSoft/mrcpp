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

/**
 *  Base class for Gaussian type functions
 *
 *  High-level overview
 *  -------------------
 *  This header declares the abstract template class Gaussian<D>, a common base
 *  for concrete Gaussian primitives used throughout MRCPP. A Gaussian here is
 *  a separable Cartesian function in D dimensions of the form
 *
 *      f(x) = coef * Π_{d=0..D-1} (x_d - pos[d])^{power[d]} * exp(-alpha[d] * (x_d - pos[d])^2),
 *
 *  where:
 *    - coef        : global scalar amplitude (double).
 *    - power[d]    : non-negative integer power of the monomial factor in dim d.
 *    - alpha[d]    : strictly positive exponent in dim d (width parameter).
 *    - pos[d]      : center coordinate in dim d.
 *
 *  The class only provides *common infrastructure* (storage, screening helpers,
 *  normalization by norm, simple algebra on prefactors, batch evaluation stubs,
 *  etc.). Concrete subclasses implement the analytic pieces that depend on the
 *  exact Gaussian flavor (e.g. GaussFunc<D>, GaussPoly<D>), such as:
 *    - evalf(...)       : the actual evaluation at a point.
 *    - evalf1D(...)     : 1D component evaluation used in tensorized loops.
 *    - calcSquareNorm() : exact L2 norm.
 *    - differentiate()  : derivative producing a polynomial × Gaussian (GaussPoly).
 *    - asGaussExp()     : expansion into sum of pure Gaussians if needed.
 *
 *  Screening and visibility
 *  ------------------------
 *  Gaussian supports optional *screening* (axis-aligned bounding boxes) to skip
 *  work on dyadic tiles that are provably negligible. See:
 *    - calcScreening(stdDeviations) : builds [A,B] bounds as ± nσ around pos.
 *    - checkScreen(n, l)            : tile-level cull test at dyadic scale n.
 *    - isVisibleAtScale(...)        : heuristic visibility vs. resolution.
 *    - isZeroOnInterval(...)        : quick interval culling via ±5σ rule.
 *
 *  Relations to other types
 *  ------------------------
 *  - GaussExp<D>: an expansion (sum) of Gaussian-like terms.
 *  - GaussFunc<D>: Gaussian with a *single* monomial factor (derived class).
 *  - GaussPoly<D>: Gaussian multiplied by a *polynomial* (derivative results).
 *
 *  Thread-safety: instances are regular value objects; no shared state.
 */

#pragma once

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <memory>

#include "MRCPP/mrcpp_declarations.h"
#include "RepresentableFunction.h"

namespace mrcpp {

/**
 * @tparam D Spatial dimension (1, 2, or 3 in MRCPP usage).
 *
 * @class Gaussian
 * @brief Abstract base for separable Cartesian Gaussians in D dimensions.
 *
 * Interface summary
 * -----------------
 *  Construction:
 *    - Gaussian(a, c, r, p)              : isotropic exponent (alpha[d]=a).
 *    - Gaussian(alpha[], c, r, p)        : anisotropic exponents per axis.
 *
 *  Core virtuals (must be implemented by derived classes):
 *    - copy()                : virtual clone (CRTP alternative).
 *    - evalf(r)              : value at point r (D-vector).
 *    - evalf1D(x, d)         : 1D factor along axis d (helper).
 *    - calcSquareNorm()      : exact ∥f∥².
 *    - asGaussExp()          : expansion into GaussExp<D>.
 *    - differentiate(dir)    : analytic derivative → GaussPoly<D>.
 *
 *  Utilities provided here:
 *    - evalf(points, values) : batch evaluation per-axis (column-wise).
 *    - calcOverlap(inp)      : ⟨this|inp⟩ via GaussExp + Obara–Saika.
 *    - periodify(period, nσ) : semi-periodic replication into a GaussExp.
 *    - calcScreening(nσ)     : build ±nσ bounds and enable screening.
 *    - checkScreen(n, l)     : dyadic tile cull test when screening on.
 *    - normalize()           : rescale by 1/∥f∥ (uses calcSquareNorm()).
 *    - multPureGauss(lhs,rhs): complete-the-square product of *pure* Gaussians
 *                              (monomials handled by derived classes).
 *    - multConstInPlace(c)   : scale the global coefficient.
 *
 *  Accessors/mutators:
 *    - get/set for coef, alpha, pos, power; toggle screen flag.
 */
template <int D> class Gaussian : public RepresentableFunction<D, double> {
public:
    /**
     * @brief Isotropic constructor.
     * @param a   Exponent value α to be replicated on all axes (α[d] = a).
     * @param c   Global scalar coefficient.
     * @param r   Center position (Coord<D>), defaults to origin.
     * @param p   Per-axis monomial powers (non-negative).
     *
     * @warning This ctor does not check positivity of @p a; callers are expected
     *          to pass α>0 (required for square integrability and σ = 1/√(2α)).
     */
    Gaussian(double a, double c, const Coord<D> &r, const std::array<int, D> &p);

    /**
     * @brief Anisotropic constructor (per-axis exponents).
     * @param a   Exponent array α[d] per axis.
     * @param c   Global scalar coefficient.
     * @param r   Center position (Coord<D>).
     * @param p   Per-axis monomial powers (non-negative).
     */
    Gaussian(const std::array<double, D> &a, double c, const Coord<D> &r, const std::array<int, D> &p);

    Gaussian<D> &operator=(const Gaussian<D> &gp) = delete; ///< Non-assignable; use clones.
    virtual Gaussian<D> *copy() const = 0;                  ///< Virtual copy (clone).
    virtual ~Gaussian() = default;

    /** @name Evaluation API (to be implemented by subclasses) */
    ///@{
    /**
     * @brief Evaluate f(r) at a D-dimensional coordinate.
     * @param r Point (Coord<D>) in physical space.
     * @return Function value f(r).
     */
    virtual double evalf(const Coord<D> &r) const = 0;

    /**
     * @brief Evaluate the *1D* separable factor along axis @p dim.
     * @param r   Coordinate along axis @p dim.
     * @param dim Axis index in [0, D-1].
     * @return g_dim(r) = (r-pos[dim])^{power[dim]} exp[-α[dim](r-pos[dim])²], possibly scaled.
     */
    virtual double evalf1D(double r, int dim) const = 0;

    /**
     * @brief Batch evaluation helper.
     * @param points Matrix (N×D): column d holds all coordinates along axis d.
     * @param values Matrix (N×D): on return, values(i,d) = evalf1D(points(i,d), d).
     *
     * @note This does *not* multiply across dimensions; it only fills the
     *       per-axis factors column-wise for later tensor products.
     */
    void evalf(const Eigen::MatrixXd &points, Eigen::MatrixXd &values) const;
    ///@}

    /** @name Integral properties and expansions */
    ///@{
    /**
     * @brief Overlap ⟨this|inp⟩ computed via GaussExp reduction and
     *        Obara–Saika 1D recurrences per dimension.
     */
    double calcOverlap(const Gaussian<D> &inp) const;

    /** @brief Exact L2 norm squared ∥f∥² (implemented by subclass). */
    virtual double calcSquareNorm() const = 0;

    /**
     * @brief Represent as a sum of Gaussians (pure or polynomial-times-Gaussian),
     *        suitable for pairwise operations; implemented by subclass.
     */
    virtual GaussExp<D> asGaussExp() const = 0;

    /**
     * @brief Create a semi-periodic expansion by replicating the function on a
     *        Cartesian lattice so that most of its mass lies within a unit cell.
     * @param period  Per-axis period lengths.
     * @param nStdDev Number of standard deviations to preserve (default 4.0).
     */
    GaussExp<D> periodify(const std::array<double, D> &period, double nStdDev = 4.0) const;
    ///@}

    /** @name Differential operators */
    ///@{
    /**
     * @brief Analytic derivative ∂/∂x_dir (Cartesian direction).
     * @param dir Axis index in [0, D-1].
     * @return A GaussPoly<D> representing the derivative (polynomial×Gaussian).
     */
    virtual GaussPoly<D> differentiate(int dir) const = 0;
    ///@}

    /** @name Screening and normalization */
    ///@{
    /**
     * @brief Build ±nσ bounds around the center on each axis and enable screening.
     *        Used to cheaply cull tiles/intervals that cannot contribute.
     * @param stdDeviations Number of standard deviations n used for the box.
     */
    void calcScreening(double stdDeviations);

    /**
     * @brief Normalize in place by dividing by the L2 norm.
     * @note Calls calcSquareNorm() from the derived class.
     */
    void normalize() {
        double norm = std::sqrt(calcSquareNorm());
        multConstInPlace(1.0 / norm);
    }
    ///@}

    /** @name Algebra on the pure Gaussian core */
    ///@{
    /**
     * @brief Complete-the-square product of two *pure* Gaussians into *this*.
     *        Polynomial factors are handled in derived types (GaussFunc/GaussPoly).
     */
    void multPureGauss(const Gaussian<D> &lhs, const Gaussian<D> &rhs);

    /** @brief Scale the global coefficient by a constant. */
    void multConstInPlace(double c) { this->coef *= c; }

    /** @brief Shorthand for multConstInPlace. */
    void operator*=(double c) { multConstInPlace(c); }
    ///@}

    /** @name Screening access */
    ///@{
    bool getScreen() const { return screen; }
    /**
     * @brief Tile-level culling test for dyadic box at scale n and translation l.
     * @return true if the box is completely outside the screening bounds and can be skipped.
     */
    bool checkScreen(int n, const int *l) const;
    ///@}

    /** @name Parameter accessors */
    ///@{
    int getPower(int i) const { return power[i]; }
    double getCoef() const { return coef; }
    double getExp(int i) const { return alpha[i]; }
    const std::array<int, D> &getPower() const { return power; }
    const std::array<double, D> &getPos() const { return pos; }
    std::array<double, D> getExp() const { return alpha; }
    ///@}

    /** @name Parameter mutators */
    ///@{
    virtual void setPow(const std::array<int, D> &power) = 0;  ///< Set all monomial powers.
    virtual void setPow(int d, int power) = 0;                  ///< Set monomial power on axis d.
    void setScreen(bool _screen) { this->screen = _screen; }    ///< Enable/disable screening flag.
    void setCoef(double cf) { this->coef = cf; }                ///< Set global coefficient.
    void setExp(double _alpha) { this->alpha.fill(_alpha); }    ///< Set isotropic exponent α[d]=_alpha.
    void setExp(const std::array<double, D> &_alpha) { this->alpha = _alpha; } ///< Set per-axis exponents.
    void setPos(const std::array<double, D> &r) { this->pos = r; }             ///< Set center coordinates.
    ///@}

    /** @brief Stream pretty-printer (delegates to virtual print()). */
    friend std::ostream &operator<<(std::ostream &o, const Gaussian<D> &gauss) { return gauss.print(o); }

    friend class GaussExp<D>; ///< Allows GaussExp to access internals when assembling expansions.

protected:
    /** @name Core parameters (POD) */
    ///@{
    bool screen;                 ///< If true, use [A,B] screening in fast checks (set via calcScreening / setScreen).
    double coef;                 /**< Global scale factor (α in the docs above). */
    std::array<int, D> power;    /**< Monomial powers per axis (non-negative integers). */
    std::array<double, D> alpha; /**< Exponents per axis (>0). Controls width: σ_d = 1/√(2 α_d). */
    Coord<D> pos;                /**< Center coordinates. */
    ///@}

    /** @name Visibility / culling helpers used by trees and projection */
    ///@{
    /**
     * @brief Heuristic visibility vs. resolution scale and quadrature sampling.
     * @param scale   Dyadic scale (tile size ~ 2^{-scale}).
     * @param nQuadPts Number of quadrature points per tile edge.
     * @return false if the Gaussian is “too narrow” to be represented at this scale.
     */
    bool isVisibleAtScale(int scale, int nQuadPts) const;

    /**
     * @brief Quick check whether the function is essentially zero on [a,b] per axis,
     *        using a ±5σ bounding rule (implementation in the .cpp).
     */
    bool isZeroOnInterval(const double *a, const double *b) const;
    ///@}

    /**
     * @brief Maximum standard deviation across axes: max_d 1/√(2 α_d).
     * @details Used by periodify() to decide how many neighboring images to include.
     */
    double getMaximumStandardDiviation() const;

    /**
     * @brief Subclass hook for stream output; should print parameters in a readable way.
     */
    virtual std::ostream &print(std::ostream &o) const = 0;
};

} // namespace mrcpp