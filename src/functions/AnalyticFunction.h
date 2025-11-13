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
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 *
 * @brief Implementation of @ref RepresentableFunction for the datatype double
 */
template <int D, typename T = double>
class AnalyticFunction : public RepresentableFunction<D, T> {
public:
    /** @brief Default constructor; leaves the callable empty */
    AnalyticFunction() = default;

    /** @brief Virtual destructor to match the base class interface */
    ~AnalyticFunction() override = default;

    /**
     * @brief Constructor with raw pointers for the bounds
     *
     * @param f   The analytic function which is evaluated in this class
     * @param a   Optional raw pointer to an array of D lower bounds (can be nullptr)
     * @param b   Optional raw pointer to an array of D upper bounds (can be nullptr)
     */
    AnalyticFunction(std::function<T(const Coord<D> &r)> f,
                     const double *a = nullptr,
                     const double *b = nullptr)
            : RepresentableFunction<D, T>(a, b)
            , func(f) {}

    /**
     * @brief Overload constructor with std::vector for the bounds
     *
     * @param f  The analytic function which is evaluated in this class
     * @param a  Vector of D lower bounds.
     * @param b  Vector of D upper bounds.
     */
    AnalyticFunction(std::function<T(const Coord<D> &r)> f,
                     const std::vector<double> &a,
                     const std::vector<double> &b)
            : AnalyticFunction(f, a.data(), b.data()) {}

    /**
     * @brief Set the analytic function to be evaluated
     * @param f New analytic function
     */
    void set(std::function<T(const Coord<D> &r)> f) { this->func = f; }

    /**
     * @brief Evaluate the analytic function at coordinate @p r.
     * @param r Coordinate where to evaluate the function
     *
     * @details Checks if the point is within bounds before evaluating
     * 
     * @return The function value at point @p r
     */
    T evalf(const Coord<D> &r) const override {
        T val = T(0);
        if (not this->outOfBounds(r)) val = this->func(r);
        return val;
    }

protected:
    std::function<T(const Coord<D> &r)> func; ///< User-provided analytic function 
};

} // namespace mrcpp