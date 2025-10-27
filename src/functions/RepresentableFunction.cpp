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
 * # RepresentableFunction (implementation)
 *
 * A lightweight base providing **optional rectangular bounds** for
 * D-dimensional functions used across MRCPP. Derived classes supply the
 * actual function evaluation; this class only manages:
 *
 * - whether a function is **bounded** or **unbounded**;
 * - storage and lifetime of lower/upper bounds `A[d]`, `B[d]`;
 * - cheap **containment tests** via @ref outOfBounds.
 *
 * ## Interval semantics
 * Bounds are interpreted as a Cartesian product of **half-open intervals**:
 *
 * \f[
 *   \prod_{d=0}^{D-1} [A_d,\; B_d)
 * \f]
 *
 * so a point is considered out of bounds if **any** coordinate is
 * `< A_d` or `>= B_d`. This convention is important for tessellations,
 * avoiding double counting on shared faces.
 *
 * ## Ownership and copying
 * - Bounds are stored in dynamically allocated arrays `A` and `B` when the
 *   function is bounded. The destructor frees them.
 * - The **copy constructor** performs a deep copy of the bounds.
 * - The **assignment operator** in this base intentionally **does not**
 *   copy bounds (a documented “no-op” that returns `*this`). If you need to
 *   copy bounds, use the copy constructor instead, or call `setBounds()`.
 *
 * Derived functors can call `outOfBounds()` prior to expensive evaluations to
 * fast-return zeros outside the active box.
 */

#include "RepresentableFunction.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @brief Construct with optional bounds.
 *
 * @param a Pointer to the lower bounds array of length D, or `nullptr` for
 *          an unbounded function.
 * @param b Pointer to the upper bounds array of length D, or `nullptr` for
 *          an unbounded function.
 *
 * If either pointer is `nullptr`, the function is marked **unbounded** and no
 * memory is allocated. Otherwise, both arrays are deep-copied and the function
 * is marked **bounded**. Each dimension is validated to satisfy `a[d] ≤ b[d]`.
 */
template <int D, typename T>
RepresentableFunction<D, T>::RepresentableFunction(const double *a, const double *b) {
    if (a == nullptr or b == nullptr) {
        this->bounded = false;
        this->A = nullptr;
        this->B = nullptr;
    } else {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
        for (int d = 0; d < D; d++) {
            if (a[d] > b[d]) { MSG_ERROR("Lower bound > Upper bound."); }
            this->A[d] = a[d];
            this->B[d] = b[d];
        }
    }
}

/**
 * @brief Copy-construct from another function, including its bounds.
 *
 * Deep-copies the bounds if @p func is bounded; otherwise keeps the new
 * function unbounded.
 */
template <int D, typename T>
RepresentableFunction<D, T>::RepresentableFunction(const RepresentableFunction<D, T> &func) {
    if (func.isBounded()) {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
        for (int d = 0; d < D; d++) {
            A[d] = func.getLowerBounds()[d];
            B[d] = func.getUpperBounds()[d];
        }
    } else {
        this->bounded = false;
        this->A = nullptr;
        this->B = nullptr;
    }
}

/**
 * @brief Assignment operator (base): **does not copy bounds**.
 *
 * This is intentionally a no-op in the base class and returns `*this`
 * unchanged. Use the copy constructor if you want an identical object
 * including bounds, or call @ref setBounds explicitly after assignment.
 *
 * @note Derived classes may extend assignment to copy additional state; the
 *       base part will still leave bounds unchanged.
 */
template <int D, typename T>
RepresentableFunction<D, T> &
RepresentableFunction<D, T>::operator=(const RepresentableFunction<D, T> &func) {
    return *this;
}

/**
 * @brief Destructor releases bound storage if allocated.
 */
template <int D, typename T>
RepresentableFunction<D, T>::~RepresentableFunction() {
    if (this->isBounded()) {
        delete[] this->A;
        delete[] this->B;
    }
    this->A = nullptr;
    this->B = nullptr;
}

/**
 * @brief Set (or overwrite) bounds.
 *
 * @param a Lower bounds array of length D (must be non-null).
 * @param b Upper bounds array of length D (must be non-null).
 *
 * - If the function was previously unbounded, storage for `A` and `B` is
 *   allocated and the function becomes bounded.
 * - Each dimension is validated to have `a[d] ≤ b[d]`.
 */
template <int D, typename T>
void RepresentableFunction<D, T>::setBounds(const double *a, const double *b) {
    if (a == nullptr or b == nullptr) { MSG_ERROR("Invalid arguments"); }
    if (not isBounded()) {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
    }
    for (int d = 0; d < D; d++) {
        if (a[d] > b[d]) { MSG_ERROR("Lower bound > Upper bound."); }
        this->A[d] = a[d];
        this->B[d] = b[d];
    }
}

/**
 * @brief Check whether a point is outside the active bounds.
 *
 * @param r D-tuple (coordinate) to test.
 * @return `true` if unambiguously out of bounds, `false` otherwise.
 *
 * Semantics: if the function is **unbounded**, this always returns `false`.
 * If bounded, it returns `true` when **any** coordinate violates the
 * half-open interval in that dimension: `r[d] < A[d]` or `r[d] >= B[d]`.
 */
template <int D, typename T>
bool RepresentableFunction<D, T>::outOfBounds(const Coord<D> &r) const {
    if (not isBounded()) { return false; }
    for (int d = 0; d < D; d++) {
        if (r[d] < getLowerBound(d)) return true;
        if (r[d] >= getUpperBound(d)) return true;
    }
    return false;
}

/* Explicit template instantiations used in MRCPP. */
template class RepresentableFunction<1, double>;
template class RepresentableFunction<2, double>;
template class RepresentableFunction<3, double>;
template class RepresentableFunction<1, ComplexDouble>;
template class RepresentableFunction<2, ComplexDouble>;
template class RepresentableFunction<3, ComplexDouble>;

} // namespace mrcpp
