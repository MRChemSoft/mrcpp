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
#include <cmath>
#include <iostream>
#include <memory>

#include "MRCPP/mrcpp_declarations.h"
#include "RepresentableFunction.h"

namespace mrcpp {

/**
 * @class Gaussian
 * 
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 *
 * @brief Represent and manipulate Gaussian-type functions (GTFs)
 * 
 * @details The Gaussian class is an abstract base class for
 * representing Gaussian-type functions (GTFs) in D dimensions.
 * A GTF is defined as
 * \f$ f(\mathbf{r}) = C \prod_{d=0}^{D-1} g_d(r_d) \f$, where
 * \f$ g_d(r_d) = (r_d - p_d)^{\alpha_d} \exp[-\beta_d (r_d - p_d)^2] \f$.
 * Here, C is a global coefficient, p_d is the center position along axis d,
 * \alpha_d is the exponent for the monomial term, and \beta_d is the exponent for the Gaussian envelope.
 * 
 * The class provides methods for evaluating the function at given points,
 * computing overlap integrals with other Gaussian functions, differentiating
 * the function, and converting to Gaussian expansions suitable for pairwise operations.
 */
template <int D> class Gaussian : public RepresentableFunction<D, double> {
public:
    /**
     * @brief Isotropic constructor (same exponent on all axes)
     * @param a   Exponent value α to be replicated on all axes (α[d] = a)
     * @param c   Global scalar coefficient
     * @param r   Center position (Coord<D>), defaults to origin
     * @param p   Per-axis monomial powers (non-negative), stored as array
     *
     * @warning This ctor does not check positivity of @p a; callers are expected
     *          to pass α>0 (required for square integrability and σ = 1/√(2α)).
     */
    Gaussian(double a, double c, const Coord<D> &r, const std::array<int, D> &p);

    /**
     * @brief Anisotropic constructor (different set of coefficients and exponents per each axis)
     * @param a   Exponent ARRAY α[d] per axis.
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
     * @brief Evaluate the gaussian f(r) at a D-dimensional coordinate
     * @param[in] r Point (Coord<D>) in physical space in the MRA box
     * @return Function value f(r).
     */
    virtual double evalf(const Coord<D> &r) const = 0;

    /**
     * @brief Evaluate the *1D* separable factor along axis @p dim
     * @param r   Coordinate along axis @p dim
     * @param dim Axis index in [0, D-1].
     * 
     * @return The value of the 1D Gaussian factor g_dim(r), dim = {0, .., D-1} -> x,y,z...
     */
    virtual double evalf1D(double r, int dim) const = 0;

    /**
     * @brief Evaluate a set of points in D dimensions, arranged in the matrix form
     * @param[in] points Matrix (N×D): column d holds all coordinates along axis d.
     * @param[out] values Matrix (N×D): on return, values(i,d) = evalf1D(points(i,d), d).
     *
     * @note This does *not* multiply across dimensions; it only fills the
     *       per-axis factors column-wise for later tensor products.
     */
    void evalf(const Eigen::MatrixXd &points, Eigen::MatrixXd &values) const;
    ///@}

    /** @name Integral properties and expansions */
    ///@{
    /**
     * @brief Overlap ⟨this|inp⟩ of two gaussians
     * 
     * @param[in] inp The other Gaussian<D> instance
     *
     * @return The value of the overlap integral ⟨this|inp⟩ as a double
     */
    double calcOverlap(const Gaussian<D> &inp) const;

    /// @return Exact L2 norm squared ∥f∥² (implemented by subclass)
    virtual double calcSquareNorm() const = 0;

    
    /// @brief Represent as a sum of Gaussians (pure or polynomial-times-Gaussian), suitable for pairwise operations; implemented by subclass
    virtual GaussExp<D> asGaussExp() const = 0;

    /** @brief Generates a GaussExp that is semi-periodic around a unit-cell
     *
     * @returns Semi-periodic version of a Gaussian around a unit-cell
     * @param[in] period: The period of the unit cell
     * @param[in] nStdDev: Number of standard diviations covered in each direction. Default 4.0
     *
     * @details nStdDev = 1, 2, 3 and 4 ensures atleast 68.27%, 95.45%, 99.73% and 99.99% of the
     * integral is conserved with respect to the integration limits.
     */
    GaussExp<D> periodify(const std::array<double, D> &period, double nStdDev = 4.0) const;
    ///@}

    /** @name Differential operators */
    ///@{
    /**
     * @brief Analytic derivative d/dx_dir (Cartesian direction) of the Gaussian
     * @param dir Axis index in [0, D-1]
     *
     * @return A GaussPoly<D> representing the derivative (polynomial×Gaussian)
     */
    virtual GaussPoly<D> differentiate(int dir) const = 0;
    ///@}

    /** @name Screening and normalization */
    ///@{
    /**
     * @brief Build ±nσ bounds around the center on each axis and enable screening
     *
     * @param stdDeviations Number of standard deviations n used for the box
     * 
     * @note Used to cheaply cull tiles/intervals that cannot contribute
     */
    void calcScreening(double stdDeviations);

    /**
     * @brief Rescale the Gaussian so that its L2 norm equals 1.
     * @note Calls calcSquareNorm() from the derived class
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

    /// @brief Scale the global coefficient by a constant
    void multConstInPlace(double c) { this->coef *= c; }

    /// @brief Shorthand for multConstInPlace
    void operator*=(double c) { multConstInPlace(c); }
    ///@}

    /** @name Screening access */
    ///@{
    bool getScreen() const { return screen; }
    /**
     * @brief Tile-level culling test for dyadic box at scale n and translation l
     * @return True if the box is completely outside the screening bounds and can be skipped
     */
    bool checkScreen(int n, const int *l) const;
    ///@}


    // some getters and setters
    /** @name Parameter accessors */
    ///@{
    int getPower(int i) const { return power[i]; }                    ///< Get monomial power on axis i
    double getCoef() const { return coef; }                   ///< Get monomial coefficient
    double getExp(int i) const { return alpha[i]; }                   ///< Get monomial exponent on axis i
    const std::array<int, D> &getPower() const { return power; }                      ///< Get monomial powers on the axis in an array
    const std::array<double, D> &getPos() const { return pos; }                   ///< Get monomial positions on the axis in an array
    std::array<double, D> getExp() const { return alpha; }                    ///< Get monomial exponent on the axis in an array
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
    bool screen;                 ///< If true, use [A,B] screening in fast checks (set via calcScreening / setScreen)
    double coef;                 ///< Global scale factor (α in the docs above)
    std::array<int, D> power;    ///< Monomial powers per axis (non-negative integers)
    std::array<double, D> alpha; ///< Exponents per axis (>0). Controls width: σ_d = 1/√(2 α_d)
    Coord<D> pos;                ///< Center coordinates
    ///@}

    /** @name Visibility / culling helpers used by trees and projection */
    ///@{
    /**
     * @brief Heuristic visibility vs. resolution scale and quadrature sampling
     * @param scale   Dyadic scale (tile size ~ 2^{-scale})
     * @param nQuadPts Number of quadrature points per tile edge
     * @return false if the Gaussian is “too narrow” to be represented at this scale
     */
    bool isVisibleAtScale(int scale, int nQuadPts) const;

    /**
     * @brief Quick check whether the function is essentially zero on [a,b] per axis,
     *        using a ±5σ bounding rule (implementation in the .cpp)
     * @param a Lower bounds array of length D
     * @param b Upper bounds array of length D
     * @return true if the function is effectively zero on [a,b]
     */
    bool isZeroOnInterval(const double *a, const double *b) const;
    ///@}

    /**
     * @brief Maximum standard deviation across axes: max_d 1/√(2 α_d).
     * @details Used by periodify() to decide how many neighboring images to include
     * 
     * @return The maximum standard deviation among all axes
     */
    double getMaximumStandardDiviation() const;

    /**
     * @brief Subclass hook for stream output; should print parameters in a readable way
     * @param o The output stream
     * @return The output stream
     */
    virtual std::ostream &print(std::ostream &o) const = 0;
};

} // namespace mrcpp