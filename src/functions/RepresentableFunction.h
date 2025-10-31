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
 * # RepresentableFunction (interface)
 *
 * Base interface for objects that can be **represented/evaluated** in the
 * multiresolution (multiwavelet) framework. Typical implementations include
 * analytic functors, Gaussian(-like) functions/expansions, polynomials and
 * function trees.
 *
 * ## Bounding box semantics
 * A function may be marked **bounded** on a Cartesian product of *half-open*
 * intervals:
 *
 *   Π_d [ A_d, B_d )
 *
 * The half-open convention prevents double counting on shared cell faces and
 * is used consistently by `outOfBounds()`. If a function is **unbounded**, its
 * bounds pointers are `nullptr` and containment checks always succeed.
 *
 * ## Lifetime & copying
 * - Bounds (arrays `A`, `B` of length `D`) are owned by the instance when set.
 * - The copy constructor **deep-copies** the bounds (if any).
 * - The assignment operator in the base class returns `*this` (does not copy
 *   bounds), leaving copying policy to derived classes if needed.
 */

#pragma once

#include <Eigen/Core>
#include <iostream>
#include <vector>

#include "MRCPP/constants.h"
#include "MRCPP/mrcpp_declarations.h"
#include "MRCPP/utils/math_utils.h"
#include "trees/NodeIndex.h"

namespace mrcpp {

/**
 * @tparam D Spatial dimension (1, 2, 3, …).
 * @tparam T Value type returned by the function (e.g. `double`,
 *           complex types, etc.).
 *
 * @brief Abstract base class for functions evaluable in the multiwavelet basis.
 *
 * The class provides **optional bounding boxes** and related helpers, while
 * deferring the actual evaluation to @ref evalf implemented by derived types.
 */
template <int D, typename T> class RepresentableFunction {
public:
    /**
     * @name Construction & assignment
     * @{
     */

    /**
     * @brief Construct with optional bounds.
     *
     * If either `a` or `b` is `nullptr`, the function is created unbounded.
     * Otherwise, `A[d]=a[d]` and `B[d]=b[d]` are deep-copied and the function
     * becomes bounded. Each dimension is validated to satisfy `a[d] ≤ b[d]`.
     *
     * @param a Lower bounds array of length `D` or `nullptr`.
     * @param b Upper bounds array of length `D` or `nullptr`.
     */
    RepresentableFunction(const double *a = nullptr, const double *b = nullptr);

    /// Convenience constructor from `std::vector` bounds.
    RepresentableFunction(const std::vector<double> &a, const std::vector<double> &b)
            : RepresentableFunction(a.data(), b.data()) {}

    /**
     * @brief Copy-construct, including bounds if present.
     *
     * Deep-copies `A` and `B` when `func` is bounded; otherwise remains unbounded.
     */
    RepresentableFunction(const RepresentableFunction<D, T> &func);

    /**
     * @brief Assignment operator (base).
     *
     * The base implementation **does not** copy bounds and simply returns `*this`.
     * Derived classes may extend this behavior to copy additional state.
     */
    RepresentableFunction<D, T> &operator=(const RepresentableFunction<D, T> &func);

    /// Virtual destructor releases bound storage if allocated.
    virtual ~RepresentableFunction();
    /** @} */

    /**
     * @brief Evaluate the function at a given point.
     * @param r Cartesian coordinate (length-`D`).
     * @returns The function value at `r`.
     *
     * Derived classes should usually check @ref outOfBounds before performing
     * expensive work and return a zero value outside the active domain.
     */
    virtual T evalf(const Coord<D> &r) const = 0;

    /**
     * @name Bounds management
     * @{
     */

    /**
     * @brief Set (or overwrite) bounds.
     *
     * Allocates and stores deep copies of `a` and `b` (length `D`) if not already
     * bounded. Validates that `a[d] ≤ b[d]` for all `d`.
     */
    void setBounds(const double *a, const double *b);

    /**
     * @brief Clear bounds and mark the function unbounded.
     *
     * After this call, @ref isBounded returns `false` and @ref outOfBounds
     * will always return `false`.
     */
    void clearBounds();

    /// @returns `true` if the function has active bounds, `false` otherwise.
    bool isBounded() const { return this->bounded; }

    /**
     * @brief Test whether a point lies outside the active bounds.
     *
     * Implements the **half-open** check for each coordinate:
     * `r[d] < A[d] || r[d] >= B[d]`. If the function is unbounded,
     * this always returns `false`.
     */
    bool outOfBounds(const Coord<D> &r) const;

    /// @returns Lower bound in dimension `d` (requires @ref isBounded).
    double getLowerBound(int d) const { return this->A[d]; }
    /// @returns Upper bound in dimension `d` (requires @ref isBounded).
    double getUpperBound(int d) const { return this->B[d]; }

    /// @returns Pointer to the lower bounds array (length `D`) or `nullptr` if unbounded.
    const double *getLowerBounds() const { return this->A; }
    /// @returns Pointer to the upper bounds array (length `D`) or `nullptr` if unbounded.
    const double *getUpperBounds() const { return this->B; }
    /** @} */

    /// @note Bridge/adapter that may require direct access to bounds.
    friend class AnalyticAdaptor<D, T>;

protected:
    /** @name Internal state
     *  @{
     */
    bool bounded; ///< `true` if the function is currently bounded.
    double *A;    ///< Lower bounds (owned; `nullptr` if unbounded).
    double *B;    ///< Upper bounds (owned; `nullptr` if unbounded).
    /** @} */

    /**
     * @brief Optional visibility hint used by some projection routines.
     * @returns `true` when the function is expected to contribute at a given scale.
     */
    virtual bool isVisibleAtScale(int scale, int nQuadPts) const { return true; }

    /**
     * @brief Optional fast zero-test on an interval (per dimension).
     * @returns `true` if the function is provably zero on `[a,b]` (component-wise).
     */
    virtual bool isZeroOnInterval(const double *a, const double *b) const { return false; }
};

/**
 * @brief Matrix-valued evaluation interface.
 *
 * A companion interface that asks an object to produce a **batch evaluation**
 * over all quadrature points associated with a tree node, returning a matrix
 * whose layout is decided by the concrete implementation.
 *
 * This is useful for high-throughput projection steps where per-point
 * overhead must be minimized.
 */
class RepresentableFunction_M {
public:
    RepresentableFunction_M() {}

    /**
     * @brief Evaluate at all points described by a node index.
     * @param nIdx Node index (scale and translation), typically defines the
     *             evaluation grid/points.
     * @returns A matrix of values (shape and semantics are implementation-defined).
     */
    virtual Eigen::MatrixXd evalf(mrcpp::NodeIndex<3> nIdx) const = 0;
};

} // namespace mrcpp