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

#include "MRCPP/constants.h"

#include "functions/Polynomial.h"

namespace mrcpp {

/**
 * @class ScalingBasis
 * @brief Abstract base class for 1D scaling-function families used in the multiwavelet basis
 *
 * @details
 * A scaling basis is a finite set of polynomials \f$ \{\varphi_k\}_{k=0}^{q-1} \f$ with \f$ q = \text{order}+1 \f$
 * that span the scaling space at level 0 for a chosen multiwavelet family.
 * Concrete families (LegendreBasis, InterpolatingBasis) derive from this class, construct and store the
 * polynomials in #funcs, populate the evaluation matrix #quadVals at quadrature nodes, and build the
 * coefficient-to-value map #cvMap and its inverse #vcMap.
 * The family is identified by #type (Legendre or Interpol as defined in constants.h) and the polynomial
 * degree by #order.
 *
 * @note Derived classes must call the base constructor, which allocates \f$ q \times q \f$ zero matrices for
 * #quadVals, #cvMap, and #vcMap, and then fill those structures before returning from their own constructor
 */
class ScalingBasis {
public:
    /**
     * @brief Construct a scaling basis descriptor
     * @param k Polynomial order (k ≥ 0)
     * @param t Family tag (e.g., Legendre or Interpol, as defined in constants.h)
     *
     * @details Stores @p k and @p t, then allocates \f$ q \times q \f$ zero matrices for #quadVals,
     * #cvMap, and #vcMap, where \f$ q = k + 1 \f$. Derived classes fill these matrices in their own
     * constructor body.
     */
    ScalingBasis(int k, int t);
    virtual ~ScalingBasis() = default;

    /**
     * @brief Evaluate all basis polynomials at an array of sample points
     * @param[in]  r    Pointer to an array of D abscissas
     * @param[out] vals Output matrix of size \f$ q \times D \f$ with
     *                  \f$ \text{vals}(k,d) = \varphi_k(r_d) \f$ for \f$ k = 0,\ldots,q-1 \f$
     *
     * @details Iterates over each column of @p vals (one per sample point) and evaluates all \f$ q \f$
     * basis polynomials stored in #funcs. Useful for projecting or evaluating on arbitrary nodes, not
     * only at the quadrature nodes.
     *
     * @note Requires @p vals.rows() == @p q
     */
    void evalf(const double *r, Eigen::MatrixXd &vals) const;

    /** @return Mutable reference to the k-th basis polynomial φ_k */
    Polynomial &getFunc(int k) { return this->funcs[k]; }
    /** @return Const reference to the k-th basis polynomial φ_k */
    const Polynomial &getFunc(int k) const { return this->funcs[k]; }

    /** @return The type of scaling basis (Legendre or Interpol; see MRCPP/constants.h) */
    int getScalingType() const { return this->type; }
    /** @return Polynomial order k */
    int getScalingOrder() const { return this->order; }
    /** @return Quadrature order q = k + 1 (one node per basis function) */
    int getQuadratureOrder() const { return this->order + 1; }

    /** @return Matrix of basis values at quadrature nodes (q × q) */
    const Eigen::MatrixXd &getQuadratureValues() const { return this->quadVals; }

    /**
     * @brief Return the coefficient-to-value or value-to-coefficient conversion map
     * @param operation Pass @c Forward (from constants.h) to obtain #cvMap (coefficients → nodal values);
     *                  any other value returns #vcMap (nodal values → coefficients)
     * @return Const reference to the selected \f$ q \times q \f$ conversion matrix
     */
    const Eigen::MatrixXd &getCVMap(int operation) const;

    /** @brief Return true if both bases have the same family type and polynomial order */
    bool operator==(const ScalingBasis &basis) const;
    /** @brief Return true if the bases differ in family type or polynomial order */
    bool operator!=(const ScalingBasis &basis) const;

    /** @brief Stream output operator; delegates to the virtual print() helper */
    friend std::ostream &operator<<(std::ostream &o, const ScalingBasis &bas) { return bas.print(o); }

protected:
    const int type;  ///< Family tag (Legendre or Interpol; see constants.h)
    const int order; ///< Polynomial order k (so \f$ q = k+1 \f$ nodes are used)

    Eigen::MatrixXd quadVals; ///< Basis values at quadrature nodes: \f$ \text{quadVals}(i,k) = \varphi_k(x_i) \f$ (size \f$ q \times q \f$)
    Eigen::MatrixXd cvMap;    ///< Coefficient-to-value map (size \f$ q \times q \f$): maps coefficient vectors to nodal values
    Eigen::MatrixXd vcMap;    ///< Value-to-coefficient map (size \f$ q \times q \f$): maps nodal values to coefficient vectors

    std::vector<Polynomial> funcs; ///< Basis polynomials \f$ \varphi_0,\ldots,\varphi_k \f$ (size \f$ q \f$)

    /**
     * @brief Polymorphic pretty-printer invoked by operator<<
     * @param[in,out] o Output stream
     * @return Reference to @p o after writing order and family name
     */
    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp