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

#include "Gaussian.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/** 
 * @class GaussFunc
 * @tparam D Spatial dimension (1, 2, or 3)
 *
 * @brief Gaussian function in D dimensions with a simple monomial in front
 *
 * - Monodimensional Gaussian (GaussFunc<1>):
 * \f$ g(x) = \alpha (x-x_0)^a e^{-\beta (x-x_0)^2} \f$
 *
 * - Multidimensional Gaussian (GaussFunc<D>):
 * \f$ G(x) = \prod_{d=1}^D g^d(x^d) \f$
 */
template <int D> class GaussFunc : public Gaussian<D> {
public:
    /** 
     * @brief Constructor which forwads to the Gaussian<D> constructor
     * @param beta      Exponent, \f$ e^{-\beta r^2} \f$
     * @param alpha     Coefficient, \f$ \alpha e^{-r^2} \f$
     * @param[in] pos   Position \f$ (x - pos[0]), (y - pos[1]), ... \f$
     * @param[in] pow   Monomial power, \f$ x^{pow[0]}, y^{pow[1]}, ... \f$
     */
    GaussFunc(double beta, double alpha, const Coord<D> &pos = {}, const std::array<int, D> &pow = {})
            : Gaussian<D>(beta, alpha, pos, pow) {}

    /** 
     * @brief Constructor which forwads to the Gaussian<D> constructor
     * @param[in] beta  List of exponents, \f$ e^{-\beta r^2} \f$
     * @param alpha     Coefficient, \f$ \alpha e^{-r^2} \f$
     * @param[in] pos   Position \f$ (x - pos[0]), (y - pos[1]), ... \f$
     * @param[in] pow   Monomial power, \f$ x^{pow[0]}, y^{pow[1]}, ... \f$
     */
    GaussFunc(const std::array<double, D> &beta,
              double alpha,
              const Coord<D> &pos = {},
              const std::array<int, D> &pow = {})
            : Gaussian<D>(beta, alpha, pos, pow) {}

    /// @brief Copy constructor.
    GaussFunc(const GaussFunc<D> &gf)
            : Gaussian<D>(gf) {}

    GaussFunc<D> &operator=(const GaussFunc<D> &rhs) = delete;

    /**
     * @brief Performs a deep copy
     * @return Pointer to a new GaussFunc<D> copy of this instance
     */
    Gaussian<D> *copy() const override;
    
    /** 
     * @brief Compute Coulomb repulsion energy between this GaussFunc and another
     * @param gf Other GaussFunc
     * @return Coulomb energy
     * 
     * @note Implemented only for D = 3
     * @note Both Gaussians must be normalized to unit charge
     * \f$ \alpha = (\beta/\pi)^{D/2} \f$ for this to be correct!
     */
    double calcCoulombEnergy(const GaussFunc<D> &rhs) const;

    /**
     * @brief Calculates the squared norm of this GaussFunc
     * @return The squared norm
     */
    double calcSquareNorm() const override;

    /**
     * @brief Evaluate the gaussian f(r) at a D-dimensional coordinate
     * @param r Point (Coord<D>) in physical space in the MRA box
     * @return Function value f(r).
     */
    double evalf(const Coord<D> &r) const override;

    /**
     * @brief Evaluate the *1D* separable factor along axis @p dim
     * @param r   Coordinate along axis @p dim
     * @param dim Axis index in [0, D-1].
     * 
     * @return The value of the 1D Gaussian factor g_dim(r), dim = {0, .., D-1} -> x,y,z...
     */
    double evalf1D(double r, int dir) const override;

    /**
     * @brief Convert this GaussFunc to a GaussExp object
     * @return A GaussExp<D> representing this GaussFunc
     */
    GaussExp<D> asGaussExp() const override;

    /**
     * @brief Analytic derivative d/dx_dir (Cartesian direction) of the GaussFunc
     * @param dir Axis index in [0, D-1]
     *
     * @return A GaussPoly<D> representing the derivative (polynomial×Gaussian)
     */
    GaussPoly<D> differentiate(int dir) const override;

    /**
     * @brief Multiplies this GaussFunc in-place with another GaussFunc
     * @param rhs The GaussFunc to multiply with
     * @note The result is stored in this GaussFunc, thus overwriting its previous values
     */
    void multInPlace(const GaussFunc<D> &rhs);

    /**
     * @brief Operator overload forwarding to multInPlace
     * @param rhs The GaussFunc to multiply with
     */
    void operator*=(const GaussFunc<D> &rhs) { multInPlace(rhs); }

    /** 
     * @brief Multiply another GaussFunc with this GaussFunc
     * @param rhs Other GaussFunc
     * @return Resulting GaussPoly<D>
     */
    GaussPoly<D> mult(const GaussFunc<D> &rhs);

    /** 
     * @brief Multiply this GaussFunc with a scalar
     * @param c Scalar to multiply
     * @returns Resulting GaussFunc<D>
     */
    GaussFunc<D> mult(double c);

    /**
     * @brief Operator overload forwarding to mult
     * @param rhs The GaussFunc to multiply with
     * @return Resulting GaussPoly<D>
     */
    GaussPoly<D> operator*(const GaussFunc<D> &rhs) { return this->mult(rhs); }

    /**
     * @brief Operator overload forwarding to mult
     * @param rhs Scalar to multiply with
     * @return Resulting GaussFunc<D>
     */
    GaussFunc<D> operator*(double c) { return this->mult(c); }

    /**
     * @brief Set the power in dimension d
     * @param d Dimension index
     * @param power Power to set
     */
    void setPow(int d, int power) override { this->power[d] = power; }

    /**
     * @brief Set the powers in all dimensions
     * @param power Array of powers to set
     */
    void setPow(const std::array<int, D> &power) override { this->power = power; }

private:
    /// @brief Print GaussFunc to output stream
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp