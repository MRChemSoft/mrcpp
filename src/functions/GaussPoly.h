/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3 of the License, or
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
#include <vector>

#include "MRCPP/mrcpp_declarations.h"

#include "Gaussian.h"
#include "Polynomial.h"

namespace mrcpp {

/** 
 * @class GaussPoly
 * @tparam D Spatial dimension (1, 2, or 3)
 * 
 * @brief Gaussian function in D dimensions with a general polynomial in front
 *
 * - Monodimensional Gaussian (GaussPoly<1>):
 *
 * \f$ g(x) = \alpha P(x-x_0) e^{-\beta (x-x_0)^2} \f$
 *
 * - Multidimensional Gaussian (GaussFunc<D>):
 *
 * \f$ G(x) = \prod_{d=1}^D g^d(x^d) \f$
 */

template <int D> class GaussPoly : public Gaussian<D> {
public:
    /** 
     * @brief Constructor
     *
     * @param beta Exponent, \f$ e^{-\beta r^2} \f$
     * @param alpha Coefficient, \f$ \alpha e^{-r^2} \f$
     * @param[in] pos Position \f$ (x - pos[0]), (y - pos[1]), ... \f$
     * @param[in] pow Max polynomial degree, \f$ P_0(x), P_1(y), ... \f$
     */
    GaussPoly(double alpha = 0.0, double coef = 1.0, const Coord<D> &pos = {}, const std::array<int, D> &power = {});

    /** 
     * @brief Constructor
     *
     * @param[in] beta List of exponents, \f$ e^{-\beta r^2} \f$
     * @param alpha Coefficient, \f$ \alpha e^{-r^2} \f$
     * @param[in] pos Position \f$ (x - pos[0]), (y - pos[1]), ... \f$
     * @param[in] pow Max polynomial degree, \f$ P_0(x), P_1(y), ... \f$
     */
    GaussPoly(const std::array<double, D> &alpha,
              double coef,
              const Coord<D> &pos = {},
              const std::array<int, D> &power = {});

    /// @brief Copy constructor.
    GaussPoly(const GaussPoly<D> &gp);

    /** 
     * @brief Construct from a GaussFunc
     * @param[in] gf: GaussFunc to convert
     */
    GaussPoly(const GaussFunc<D> &gf);

    GaussPoly<D> &operator=(const GaussPoly<D> &gp) = delete;

    /**
     * @brief Performs a deep copy
     * @return Pointer to a new GaussFunc<D> copy of this instance
     */
    Gaussian<D> *copy() const override;

    ~GaussPoly();

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
    double evalf1D(double r, int dim) const override;

    /**
     * @brief Convert this GaussPoly to a GaussExp object
     * @return A GaussExp<D> representing this GaussPoly
     */
    GaussExp<D> asGaussExp() const override;

    /// @warning This method is currently not implemented.
    GaussPoly differentiate(int dir) const override;

    /// @warning This method is currently not implemented.
    void multInPlace(const GaussPoly<D> &rhs);

    /** @brief In-place product operator (delegates to @ref multInPlace). */
    void operator*=(const GaussPoly<D> &rhs) { multInPlace(rhs); }

    /// @warning This method is currently not implemented.
    GaussPoly<D> mult(const GaussPoly<D> &rhs);

    /** 
     * @brief Multiply this GaussPoly with a scalar
     * @param c Scalar to multiply
     * @returns Resulting GaussPoly<D>
     */
    GaussPoly<D> mult(double c);

    /**
     * @brief Operator overload forwarding to mult
     * @param rhs The GaussPoly to multiply with
     * @return Resulting GaussPoly<D>
     * @warning @ref mult is currently not implemented.
     */
    GaussPoly<D> operator*(const GaussPoly<D> &rhs) { return mult(rhs); }

    /**
     * @brief Operator overload forwarding to mult
     * @param rhs Scalar to multiply with
     * @return Resulting GaussPoly<D>
     */
    GaussPoly<D> operator*(double c) { return mult(c); }

    /**
     * @brief Returns the polynomial coefficients in a specified dimension
     * @param i Dimension index
     * @return The Eigen vector of coefficients
     */
    const Eigen::VectorXd &getPolyCoefs(int i) const { return poly[i]->getCoefs(); }

    /** 
     * @brief Returns the polynomial coefficients in a specified dimension
     * @param i Dimension index
     * @return The Eigen vector of coefficients
     */
    Eigen::VectorXd &getPolyCoefs(int i) { return poly[i]->getCoefs(); }

    /**
     * @brief Returns the Polynomial in a specified dimension
     * @param i Dimension index
     * @return The Polynomial reference
     */
    const Polynomial &getPoly(int i) const { return *poly[i]; }

    /** 
     * @brief Returns the Polynomial in a specified dimension
     * @param i Dimension index
     * @return The Polynomial reference
     */
    Polynomial &getPoly(int i) { return *poly[i]; }

    /**
     * @brief Set the power in dimension d
     * @param d Dimension index
     * @param power Power to set
     */
    void setPow(int d, int pow) override;

    /**
     * @brief Set the powers in all dimensions
     * @param power Array of powers to set
     */
    void setPow(const std::array<int, D> &pow) override;

    /** 
     * @brief Set polynomial in given dimension
     * @param d Cartesian direction
     * @param[in] poly Polynomial to set
     */
    void setPoly(int d, Polynomial &poly);

private:
    Polynomial *poly[D]; ///< Per-axis polynomial factors

    /**
     * @brief Recursive helper function to fill coefficient and power vectors for all terms
     * @param[out] coefs Vector to fill with coefficients
     * @param[out] power Vector to fill with power arrays
     * @param pow Current power array being built
     * @param dir Current dimension being processed
     */
    void fillCoefPowVector(std::vector<double> &coefs, std::vector<int *> &power, int pow[D], int dir) const;

    /**
     * @brief Recursive helper function to fill coefficient and power vectors for all terms
     * @param[out] coefs Vector to fill with coefficients
     * @param[out] power Vector to fill with power arrays
     * @param pow Current power array being built
     * @param dir Current dimension being processed
     */
    void fillCoefPowVector(std::vector<double> &coefs,
                           std::vector<int *> &power,
                           std::array<int, D> &pow,
                           int dir) const;

    /// @brief Print GaussFunc to output stream
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp