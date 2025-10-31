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

#include <functional>
#include <vector>

#include "RepresentableFunction.h"

namespace mrcpp {

/**
 * @class AnalyticFunction
 * @tparam D Spatial dimension (1, 2, 3, …).
 * @tparam T Numeric value type (defaults to double).
 *
 * @brief Thin adapter that turns a C++ callable `std::function<T(const Coord<D>&)>`
 *        into a @ref RepresentableFunction suitable for MRCPP algorithms.
 *
 * Motivation
 * ----------
 * Many MRCPP routines operate on the abstract interface `RepresentableFunction<D,T>`
 * (which provides domain bounds and an `evalf()` method). `AnalyticFunction` lets
 * users plug in any analytic lambda or function pointer without writing a full
 * derived class.
 *
 * Domain handling
 * ---------------
 * The base class @ref RepresentableFunction stores lower/upper bounds for each
 * coordinate dimension. `AnalyticFunction::evalf` first checks
 * `RepresentableFunction::outOfBounds(r)` and **returns 0** for points outside
 * the domain; otherwise it forwards to the user-supplied callable.
 *
 * Typical usage
 * -------------
 * @code
 * using F = AnalyticFunction<2>;
 * std::vector<double> a = {0.0, 0.0};
 * std::vector<double> b = {1.0, 2.0};
 * F f(
 *   [](const Coord<2>& x) -> double {
 *       // x[0] = x, x[1] = y
 *       return std::sin(x[0]) * std::exp(-x[1]);
 *   },
 *   a, b
 * );
 * Coord<2> p; p[0] = 0.3; p[1] = 1.5;
 * double v = f.evalf(p); // evaluates lambda if p within [a,b]
 * @endcode
 *
 * Thread-safety
 * -------------
 * `AnalyticFunction` itself holds only an immutable std::function after construction.
 * It is safe to call `evalf` concurrently *iff your callable is thread-safe* and
 * does not mutate shared state.
 */
template <int D, typename T = double>
class AnalyticFunction : public RepresentableFunction<D, T> {
public:
    /** @brief Default constructor; leaves the callable empty. */
    AnalyticFunction() = default;

    /** @brief Virtual destructor to match the base class interface. */
    ~AnalyticFunction() override = default;

    /**
     * @brief Construct with a callable and optional raw-pointer bounds.
     *
     * @param f   Callable of signature `T(const Coord<D>&)`.
     * @param a   Optional pointer to an array of D lower bounds (can be nullptr).
     * @param b   Optional pointer to an array of D upper bounds (can be nullptr).
     *
     * The bounds are forwarded to the @ref RepresentableFunction base; if both
     * are nullptr the base uses its defaults (implementation-defined).
     */
    AnalyticFunction(std::function<T(const Coord<D> &r)> f,
                     const double *a = nullptr,
                     const double *b = nullptr)
            : RepresentableFunction<D, T>(a, b)
            , func(f) {}

    /**
     * @brief Construct with a callable and STL vector bounds.
     *
     * @param f Callable of signature `T(const Coord<D>&)`.
     * @param a Vector of D lower bounds.
     * @param b Vector of D upper bounds.
     *
     * Convenience overload that forwards raw pointers of the vectors to the
     * other constructor. The vectors must have length D.
     */
    AnalyticFunction(std::function<T(const Coord<D> &r)> f,
                     const std::vector<double> &a,
                     const std::vector<double> &b)
            : AnalyticFunction(f, a.data(), b.data()) {}

    /**
     * @brief Replace the underlying callable at runtime.
     *
     * @param f New callable `T(const Coord<D>&)`.
     *
     * No synchronization is performed; if other threads may call `evalf`
     * concurrently, arrange external synchronization.
     */
    void set(std::function<T(const Coord<D> &r)> f) { this->func = f; }

    /**
     * @brief Evaluate the function at coordinate @p r.
     *
     * Behavior:
     *  - If @p r lies outside the domain bounds (per `outOfBounds(r)`), return 0.
     *  - Otherwise, return `func(r)`.
     *
     * @note Returning 0 outside the domain is consistent with how many MRCPP
     *       integrators and projectors treat functions on bounded supports.
     */
    T evalf(const Coord<D> &r) const override {
        T val = T(0);
        if (not this->outOfBounds(r)) val = this->func(r);
        return val;
    }

protected:
    /**
     * @brief Stored analytic callable.
     *
     * The signature uses `Coord<D>` (MRCPP’s fixed-size coordinate array).
     * The callable should be side-effect free or externally synchronized if
     * used from multiple threads.
     */
    std::function<T(const Coord<D> &r)> func;
};

} // namespace mrcpp