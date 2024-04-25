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

#include <Eigen/Core>

#include "Gaussian.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/** @class GaussFunc
 *
 * @brief Gaussian function in D dimensions with a simple monomial in front
 *
 * - Monodimensional Gaussian (GaussFunc<1>):
 *
 * \f$ g(x) = \alpha (x-x_0)^a e^{-\beta (x-x_0)^2} \f$
 *
 * - Multidimensional Gaussian (GaussFunc<D, T>):
 *
 * \f$ G(x) = \prod_{d=1}^D g^d(x^d) \f$
 */

template <int D, typename T> class GaussFunc : public Gaussian<D, T> {
public:
    /** @returns New GaussFunc object
     *  @param[in] beta: Exponent, \f$ e^{-\beta r^2} \f$
     *  @param[in] alpha: Coefficient, \f$ \alpha e^{-r^2} \f$
     *  @param[in] pos: Position \f$ (x - pos[0]), (y - pos[1]), ... \f$
     *  @param[in] pow: Monomial power, \f$ x^{pow[0]}, y^{pow[1]}, ... \f$
     */
    GaussFunc(double beta, double alpha, const Coord<D> &pos = {}, const std::array<int, D> &pow = {})
            : Gaussian<D, T>(beta, alpha, pos, pow) {}
    GaussFunc(const std::array<double, D> &beta,
              double alpha,
              const Coord<D> &pos = {},
              const std::array<int, D> &pow = {})
            : Gaussian<D, T>(beta, alpha, pos, pow) {}
    GaussFunc(const GaussFunc<D, T> &gf)
            : Gaussian<D, T>(gf) {}
    GaussFunc<D, T> &operator=(const GaussFunc<D, T> &rhs) = delete;
    Gaussian<D, T> *copy() const override;

    double calcCoulombEnergy(const GaussFunc<D, T> &rhs) const;
    double calcSquareNorm() const override;

    T evalf(const Coord<D> &r) const override;
    T evalf1D(double r, int dir) const override;

    GaussExp<D, T> asGaussExp() const override;
    GaussPoly<D, T> differentiate(int dir) const override;

    void multInPlace(const GaussFunc<D, T> &rhs);
    void operator*=(const GaussFunc<D, T> &rhs) { multInPlace(rhs); }
    GaussPoly<D, T> mult(const GaussFunc<D, T> &rhs);
    GaussFunc<D, T> mult(double c);
    GaussPoly<D, T> operator*(const GaussFunc<D, T> &rhs) { return this->mult(rhs); }
    GaussFunc<D, T> operator*(double c) { return this->mult(c); }

    void setPow(int d, int power) override { this->power[d] = power; }
    void setPow(const std::array<int, D> &power) override { this->power = power; }

private:
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp
