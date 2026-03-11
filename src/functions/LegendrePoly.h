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

#include "Polynomial.h"

namespace mrcpp {

 /**
  * @class LegendrePoly
  * @brief Class defining a Legendre polynomial of degree k
  */
class LegendrePoly final : public Polynomial {
public:
    /**
     * @brief Construct and compute a Legendre polynomial of degree k
     * @param k Degree (order) of the Legendre polynomial 
     * @param n Dilation factor (applied after translation)
     * @param l Translation (applied before dilation)
     */
    LegendrePoly(int k, double n = 1.0, double l = 0.0);

    /**
     * @brief Evaluate value and first derivative of this Legendre polynomial in x
     * @param x External evaluation point
     * @return Value and first derivative as an Eigen::Vector2d
     */
    Eigen::Vector2d firstDerivative(double x) const;

    /**
     * @brief Evaluate second derivative of this Legendre polynomial in x
     * @param x External evaluation point
     * @return Value, first and second derivative as an Eigen::Vector3d
     */
    Eigen::Vector3d secondDerivative(double x) const;

private:
    /** 
     * @brief Recursively compute the Legendre polynomial of order k on interval [-1,1]
     * @param k Order of the Legendre polynomial
     */
    void computeLegendrePolynomial(int k);
};

} // namespace mrcpp