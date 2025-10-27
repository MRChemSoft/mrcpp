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
 * ScalingBasis provides common functionality shared by concrete scaling bases
 * (e.g., LegendreBasis, InterpolatingBasis). It stores:
 *   • the basis family tag (type) and polynomial order (order),
 *   • the list of basis polynomials (this->funcs) managed in the headers,
 *   • matrices used to convert between coefficient and nodal-value spaces:
 *       - quadVals : basis evaluated at quadrature nodes (q×q),
 *       - cvMap    : coefficient → value map at nodes (q×q),
 *       - vcMap    : value → coefficient map at nodes (q×q).
 *
 * Responsibilities in this file:
 *   - Construct and size the conversion matrices based on the quadrature order.
 *   - Provide a generic evaluator `evalf` to sample the basis at arbitrary points.
 *   - Expose the proper conversion map (cvMap or vcMap) given an operation flag.
 *   - Define equality/inequality operators and a simple printer.
 *
 * Conventions:
 *   - q := getQuadratureOrder() is the number of basis functions and nodes.
 *   - `type` identifies the scaling family (Legendre vs Interpol); codes live
 *     in shared headers.
 *   - `Forward` indicates coefficient→value mapping; anything else selects the
 *     reverse map (value→coefficient).
 */

#include "ScalingBasis.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @brief Construct a base scaling space for family @p t and order @p k.
 *
 * Initializes:
 *   - type/order (with a minimal validity check on order),
 *   - square q×q matrices quadVals, cvMap, vcMap filled with zeros,
 *     where q = getQuadratureOrder() is determined by the concrete basis.
 *
 * Concrete derived classes are responsible for:
 *   - populating `funcs` with q basis polynomials,
 *   - filling `quadVals`,
 *   - building `cvMap` and `vcMap`.
 */
ScalingBasis::ScalingBasis(int k, int t)
        : type(t)
        , order(k) {
    if (this->order < 0) MSG_ABORT("Invalid scaling order");
    int q_order = getQuadratureOrder();
    this->quadVals = Eigen::MatrixXd::Zero(q_order, q_order); // basis@nodes
    this->cvMap    = Eigen::MatrixXd::Zero(q_order, q_order); // coeff → values
    this->vcMap    = Eigen::MatrixXd::Zero(q_order, q_order); // values → coeff
}

/**
 * @brief Evaluate each basis polynomial at a set of points.
 *
 * @param[in]  r    Pointer to an array of length D containing evaluation points.
 * @param[out] vals Matrix of size (K × D) where:
 *                    - K must equal the number of basis functions (funcs.size()),
 *                    - column d receives the vector [ φ_0(r[d]), …, φ_{K-1}(r[d]) ]^T.
 *
 * Precondition:
 *   - vals.rows() == funcs.size(). If not, an error is reported.
 *
 * Notes:
 *   - The layout is "basis index in rows, sample index in columns".
 *   - getFunc(k) returns the k-th polynomial; evalf(x) evaluates it at x.
 */
void ScalingBasis::evalf(const double *r, Eigen::MatrixXd &vals) const {
    if (vals.rows() != this->funcs.size()) MSG_ERROR("Invalid argument");

    for (int d = 0; d < vals.cols(); d++) {
        for (int k = 0; k < vals.rows(); k++) {
            vals(k, d) = getFunc(k).evalf(r[d]);
        }
    }
}

/**
 * @brief Retrieve the appropriate coefficient/value conversion map.
 *
 * @param operation If equal to Forward, return cvMap (coeff → values at nodes).
 *                  Otherwise, return vcMap (values at nodes → coeff).
 *
 * The precise enum/integer value of Forward is defined in shared headers.
 * Derived classes ensure cvMap/vcMap are properly populated in their init code.
 */
const Eigen::MatrixXd &ScalingBasis::getCVMap(int operation) const {
    if (operation == Forward) {
        return this->cvMap;
    } else {
        return this->vcMap;
    }
}

/**
 * @brief Two scaling bases are equal iff both family type and order match.
 */
bool ScalingBasis::operator==(const ScalingBasis &basis) const {
    if (this->type != basis.type) return false;
    if (this->order != basis.order) return false;
    return true;
}

/**
 * @brief Negation of operator== (true if type or order differs).
 */
bool ScalingBasis::operator!=(const ScalingBasis &basis) const {
    if (this->type != basis.type) return true;
    if (this->order != basis.order) return true;
    return false;
}

/**
 * @brief Stream printer with a minimal summary (order and family name).
 *
 * Prints:
 *   - "polynomial order      : <order>"
 *   - "polynomial type       : <Legendre|Interpolating|Unknown>"
 */
std::ostream &ScalingBasis::print(std::ostream &o) const {
    o << " polynomial order      : " << getScalingOrder() << std::endl;
    if (getScalingType() == Legendre) {
        o << " polynomial type       : Legendre";
    } else if (getScalingType() == Interpol) {
        o << " polynomial type       : Interpolating";
    } else {
        o << " polynomial type       : Unknown";
    }
    return o;
}

} // namespace mrcpp