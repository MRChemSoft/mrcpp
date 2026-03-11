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
 * @brief Abstract base for scaling-function families (Legendre, Interpolating).
 *
 * What this class represents
 * --------------------------
 * A *scaling basis* is a finite set of 1D polynomials {φ_k}_{k=0..order}
 * that span the scaling space at level 0 for a given multiwavelet family.
 * Concrete families (e.g., LegendreBasis, InterpolatingBasis) derive from
 * this class and:
 *   • construct and store the polynomials in `funcs`,
 *   • populate the evaluation matrix at quadrature nodes `quadVals`,
 *   • build coefficient↔value conversion maps `cvMap` / `vcMap`.
 *
 * Dimensions and conventions
 * --------------------------
 * - order := polynomial degree cutoff (≥ 0).
 * - Quadrature order q = order + 1 (one node per basis function).
 * - `quadVals` is q×q with layout: rows = nodes, cols = basis index.
 * - `cvMap`  maps coefficient vectors → nodal values (Forward).
 * - `vcMap`  maps nodal values → coefficient vectors (Backward).
 *
 * Responsibilities provided here
 * ------------------------------
 * - Store family `type` (Legendre or Interpol, defined in constants.h) and `order`.
 * - Provide access to basis polynomials and to the conversion matrices.
 * - Offer a generic evaluator to sample the basis at arbitrary points.
 * - Define equality operators (same family and order).
 *
 * Notes for implementers of derived classes
 * -----------------------------------------
 * - Call the base ctor with (k, t). It sizes `quadVals`, `cvMap`, `vcMap`
 *   to q×q zeros; you must fill them in your implementation (.cpp).
 * - Push back exactly q polynomials into `funcs` in the order k = 0..order.
 */
class ScalingBasis {
public:
    /**
     * @brief Construct a base scaling space descriptor.
     * @param k Polynomial order (k ≥ 0).
     * @param t Family tag (e.g., Legendre or Interpol).
     *
     * Effects (implemented in the .cpp):
     *  - Stores @p t, @p k.
     *  - Allocates q×q zero matrices for `quadVals`, `cvMap`, and `vcMap`,
     *    where q = k + 1.
     *  - Derived classes then fill these structures.
     */
    ScalingBasis(int k, int t);
    virtual ~ScalingBasis() = default;

    /**
     * @brief Evaluate all basis polynomials at D sample points.
     * @param r    Pointer to array of D abscissas.
     * @param vals Output matrix of size (q × D) with
     *             vals(k, d) = φ_k( r[d] ), k = 0..q-1.
     *
     * Precondition:
     *  - vals.rows() == funcs.size() == q.
     *
     * Remarks:
     *  - Column-major Eigen storage is irrelevant here; we just fill entries.
     *  - Useful for projecting/evaluating on arbitrary nodes (not only quadrature).
     */
    void evalf(const double *r, Eigen::MatrixXd &vals) const;

    /** @return Mutable reference to the k-th basis polynomial φ_k. */
    Polynomial &getFunc(int k) { return this->funcs[k]; }
    /** @return Const reference to the k-th basis polynomial φ_k. */
    const Polynomial &getFunc(int k) const { return this->funcs[k]; }

    /** @return The type of scaling basis (Legendre or Interpol; see MRCPP/constants.h) */
    int getScalingType() const { return this->type; }
    /** @return Polynomial order k. */
    int getScalingOrder() const { return this->order; }
    /** @return Quadrature order q = k + 1 (one node per basis function). */
    int getQuadratureOrder() const { return this->order + 1; }

    /** @return Matrix of basis values at quadrature nodes (q × q). */
    const Eigen::MatrixXd &getQuadratureValues() const { return this->quadVals; }

    /**
     * @brief Access the coefficient/value conversion map.
     * @param operation Use `Forward` (from constants.h) for coeff→value,
     *                  anything else selects value→coeff.
     * @return const reference to `cvMap` (Forward) or `vcMap` (Backward).
     */
    const Eigen::MatrixXd &getCVMap(int operation) const;

    /** @brief Equality iff same family type and polynomial order. */
    bool operator==(const ScalingBasis &basis) const;
    /** @brief Inequality iff family type or polynomial order differs. */
    bool operator!=(const ScalingBasis &basis) const;

    /**
     * @brief Stream print helper (delegates to virtual print()).
     * Prints order and a human-readable family name.
     */
    friend std::ostream &operator<<(std::ostream &o, const ScalingBasis &bas) { return bas.print(o); }

protected:
    /** @brief Family tag (Legendre or Interpol). */
    const int type;
    /** @brief Polynomial order k. */
    const int order;

    /** @brief Basis values at quadrature points: quadVals(i,k) = φ_k(x_i). */
    Eigen::MatrixXd quadVals; // function values at quadrature pts

    /** @brief Coefficient → value (at nodes) linear map (q × q). */
    Eigen::MatrixXd cvMap;    // coef-value transformation matrix

    /** @brief Value (at nodes) → coefficient linear map (q × q). */
    Eigen::MatrixXd vcMap;    // value-coef transformation matrix

    /** @brief List of basis polynomials φ_0..φ_k (size q). */
    std::vector<Polynomial> funcs;

    /**
     * @brief Polymorphic pretty-printer called by operator<<.
     * Concrete bases may override to append family-specific info.
     */
    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp