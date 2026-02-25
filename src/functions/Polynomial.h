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

#include <vector>

#include <Eigen/Core>

#include "RepresentableFunction.h"

namespace mrcpp {

/**
 * @class Polynomial
 * 
 * @brief Base class for general polynomials
 * 
 * @details The Polynomial class(es) are not implemented in the
 * most efficient manner, because they are only evaluated a fixed
 * number of times in a few predefined points, and all other
 * evaluations are done by linear transformations. PolynomialCache
 * implements the fast, and static const versions of the various
 * 4Polynomials.
 */
class Polynomial : public RepresentableFunction<1, double> {
public:
    /**
     * @brief Construct polynomial of order zero with given bounds
     * @param k Order of the polynomial
     * @param a Lower bound in x as raw pointer
     * @param b Upper bound in x as raw pointer
     */
    Polynomial(int k = 0, const double *a = nullptr, const double *b = nullptr);

    /**
     * @brief Construct polynomial of order k with given bounds
     * @param k Order of the polynomial
     * @param a Lower bound in x as vector
     * @param b Upper bound in x as vector
     */
    Polynomial(int k, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(k, a.data(), b.data()) {}

    /**
     * @brief Construct polynomial with given coefficient, order and bounds
     * @param c Coefficient of the polynomial
     * @param k Order of the polynomial
     * @param a Lower bound in x as raw pointer
     * @param b Upper bound in x as raw pointer
     */
    Polynomial(double c, int k = 0, const double *a = nullptr, const double *b = nullptr);

    /**
     * @brief Construct polynomial with given coefficient, order and bounds
     * @param c Coefficient of the polynomial
     * @param k Order of the polynomial
     * @param a Lower bound in x as vector
     * @param b Upper bound in x as vector
     */
    Polynomial(double c, int k, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(c, k, a.data(), b.data()) {}

    /**
     * @brief Construct polynomial with given coefficient vector and bounds
     * @param c Coefficient vector
     * @param a Lower bound in x as raw pointer
     * @param b Upper bound in x as raw pointer
     */
    Polynomial(const Eigen::VectorXd &c, const double *a = nullptr, const double *b = nullptr);

    /**
     * @brief Construct polynomial with given coefficient vector and bounds
     * @param c Coefficient vector
     * @param a Lower bound in x as vector
     * @param b Upper bound in x as vector
     */
    Polynomial(const Eigen::VectorXd &c, const std::vector<double> &a, const std::vector<double> &b)
            : Polynomial(c, a.data(), b.data()) {}

    /** @brief Copy constructor */    
    Polynomial(const Polynomial &poly);
    /** @brief Assignment operator, copies oly the function, not its bounds */
    Polynomial &operator=(const Polynomial &poly);
    /** @brief Virtual destructor */
    virtual ~Polynomial() = default;

    /**
     * @brief Evaluate scaled and translated polynomial
     * @param x External evaluation point
     * @return The polynomial value at x
     */
    double evalf(double x) const;
    
    /**
     * @brief Evaluate scaled and translated polynomial at a given point
     * @param r 1D-Cartesian coordinate
     * @return The polynomial value at r
     */
    double evalf(const Coord<1> &r) const { return evalf(r[0]); }

    double getScaledLowerBound() const; ///< @return The actual scaled lower bound
    double getScaledUpperBound() const; ///< @return The actual scaled upper bound

    void normalize(); ///< @brief Divide by norm of (bounded) polynomial

    /**
     * @brief Calculated squared L2 norm of the (bounded) polynomial
     * @return Squared L2 norm, -1 if unbounded
     */
    double calcSquareNorm();

    double getTranslation() const { return this->L; }           ///< @return Current translation
    double getDilation()   const { return this->N; }            ///< @return Current dilation

    void setDilation(double n)    { this->N = n; }              ///< @brief Set dilation factor N
    void setTranslation(double l) { this->L = l; }              ///< @brief Set translation L
    void dilate(double n)         { this->N *= n; }             ///< @brief Dilate by factor n
    void translate(double l)      { this->L += this->N * l; }   ///< @brief Translate by l

    int size() const { return this->coefs.size(); }             ///< @return The size of the coefficient vector
    int getOrder() const;                                       ///< @return The order of the highest non-zero coefficient
    
    void clearCoefs() { this->coefs = Eigen::VectorXd::Zero(1); }               ///< @brief Clear all coefficients
    void setZero() { this->coefs = Eigen::VectorXd::Zero(this->coefs.size()); } ///< @brief Set all coefficients to zero
    void setCoefs(const Eigen::VectorXd &c) { this->coefs = c; }                ///< @brief Replace the coefficient vector with a new one

    Eigen::VectorXd &getCoefs() { return this->coefs; }                 ///< @return The coefficient vector                 
    const Eigen::VectorXd &getCoefs() const { return this->coefs; }     ///< @return The coefficient vector (const version)

    /**
     * @brief Calculates the derivative \f$ Q = dP/dx \f$ of this polynomial
     * @return The derivative polynomial Q
     */
    Polynomial calcDerivative() const;

    /**
     * @brief Calculates the indefinite integral \f$ Q = \int P\,dx \f$ of this polynomial, with constant = 0
     * @return The indefinite integral polynomial Q
     */
    Polynomial calcAntiDerivative() const;

    /**
     * @brief Calculates the derivative \f$ P \leftarrow dP/dx \f$ of this polynomial in-place
     * @details Replaces the current polynomial with its derivative, i.e. \f$ P \leftarrow dP/dx \f$.
     */
    void calcDerivativeInPlace();
    
    /**
     * @brief Calculates the indefinite integral \f$ P \leftarrow \int P\,dx \f$ of this polynomial in-place, with constant = 0
     * @details Replaces the current polynomial with its indefinite integral, i.e. \f$ P \leftarrow \int P\,dx \f$, with integration constant set to zero.
     */
    void calcAntiDerivativeInPlace();

    /**
     * @brief Calculates the analytical integral of P on [a, b]
     * @param a Lower bound of the integration interval, defaults to the polynomial's lower bound
     * @param b Upper bound of the integration interval, defaults to the polynomial's upper bound
     * @return The integral of P on [a, b]
     */
    double integrate(const double *a = 0, const double *b = 0) const;
 
    /**
     * @brief Analytically calculates the inner product of this polynomial with another one
     * @param Q The other polynomial
     * @return The inner product of the two polynomials
     */
    double innerProduct(const Polynomial &Q) const;

    /**
     * @brief In-place sum \f$ P \leftarrow P + c\,Q \f$.
     * @param c Scalar multiplier for Q
     * @param Q The polynomial to be added to P
     */
    void addInPlace(double c, const Polynomial &Q);

    /**
     * @brief Sum \f$ R = P + c\,Q \f$.
     * @param c Scalar multiplier for Q
     * @param Q The polynomial to be added to P
     * @return The resulting polynomial
     */
    Polynomial add(double c, const Polynomial &Q) const;

    /** 
     * @brief Scalar product of Polynomial with c
     * @param c The scalar multiplier
     * @return The resulting polynomial
     */
    Polynomial operator*(double c) const;

    /**
     * @brief Product of two Polynomials
     * @param Q The other polynomial
     * @return The resulting (unbounded) polynomial
     */
    Polynomial operator*(const Polynomial &Q) const;

    /**
     * @brief Sum two Polynomials
     * @param Q The other polynomial
     * @return The resulting polynomial
     */
    Polynomial operator+(const Polynomial &Q) const { return add(1.0, Q); }
    
    /**
     * @brief Difference of two Polynomials
     * @param Q The other polynomial
     * @return The resulting polynomial
     */
    Polynomial operator-(const Polynomial &Q) const { return add(-1.0, Q); }

    /**
     * @brief In-place scalar product.
     * @param c The scalar multiplier
     * @return Reference to the modified polynomial
     */
    Polynomial &operator*=(double c);
    
    /**
     * @brief In-place product of two Polynomials
     * @param Q The other polynomial
     * @return Reference to the modified polynomial
     */
    Polynomial &operator*=(const Polynomial &Q);

    /**
     * @brief In-place sum of two Polynomials
     * @param Q The other polynomial
     * @return Reference to the modified polynomial
     */
    Polynomial &operator+=(const Polynomial &Q);

    /**
     * @brief In-place difference of two Polynomials
     * @param Q The other polynomial
     * @return Reference to the modified polynomial
     */
    Polynomial &operator-=(const Polynomial &Q);

protected:
    double N;              ///< Dilation coefficient
    double L;              ///< Translation coefficient
    Eigen::VectorXd coefs; ///< Expansion coefficients
};

} // namespace mrcpp