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
 * - Multidimensional Gaussian (GaussFunc<D>):
 *
 * \f$ G(x) = \prod_{d=1}^D g^d(x^d) \f$
 */

template <int D> class GaussFunc : public Gaussian<D> {
public:
    /** @returns New GaussFunc object
     *  @param[in] beta: Exponent, \f$ e^{-\beta r^2} \f$
     *  @param[in] alpha: Coefficient, \f$ \alpha e^{-r^2} \f$
     *  @param[in] pos: Position \f$ (x - pos[0]), (y - pos[1]), ... \f$
     *  @param[in] pow: Monomial power, \f$ x^{pow[0]}, y^{pow[1]}, ... \f$
     */
    GaussFunc(double beta, double alpha, const Coord<D> &pos = {}, const std::array<int, D> &pow = {})
            : Gaussian<D>(beta, alpha, pos, pow) {}
    GaussFunc(const std::array<double, D> &beta,
              double alpha,
              const Coord<D> &pos = {},
              const std::array<int, D> &pow = {})
            : Gaussian<D>(beta, alpha, pos, pow) {}
    GaussFunc(const GaussFunc<D> &gf)
            : Gaussian<D>(gf) {}
    GaussFunc<D> &operator=(const GaussFunc<D> &rhs) = delete;
    Gaussian<D> *copy() const override;

    double calcCoulombEnergy(const GaussFunc<D> &rhs) const;
    double calcSquareNorm() const override;

    double evalf(const Coord<D> &r) const override;
    double evalf1D(double r, int dir) const override;

    GaussExp<D> asGaussExp() const override;
    GaussPoly<D> differentiate(int dir) const override;

    void multInPlace(const GaussFunc<D> &rhs);
    void operator*=(const GaussFunc<D> &rhs) { multInPlace(rhs); }
    GaussPoly<D> mult(const GaussFunc<D> &rhs);
    GaussFunc<D> mult(double c);
    GaussPoly<D> operator*(const GaussFunc<D> &rhs) { return this->mult(rhs); }
    GaussFunc<D> operator*(double c) { return this->mult(c); }

    void setPow(int d, int power) override { this->power[d] = power; }
    void setPow(const std::array<int, D> &power) override { this->power = power; }

private:
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp
