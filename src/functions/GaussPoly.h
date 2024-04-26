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
#include <vector>

#include "MRCPP/mrcpp_declarations.h"

#include "Gaussian.h"
#include "Polynomial.h"

namespace mrcpp {

/** @class GaussPoly
 *
 * @brief Gaussian function in D dimensions with a general polynomial in front
 *
 * - Monodimensional Gaussian (GaussPoly<1>):
 *
 * \f$ g(x) = \alpha P(x-x_0) e^{-\beta (x-x_0)^2} \f$
 *
 * - Multidimensional Gaussian (GaussFunc<D, T>):
 *
 * \f$ G(x) = \prod_{d=1}^D g^d(x^d) \f$
 */

template <int D, typename T> class GaussPoly : public Gaussian<D, T> {
public:
    GaussPoly(double alpha = 0.0, double coef = 1.0, const Coord<D> &pos = {}, const std::array<int, D> &power = {});
    GaussPoly(const std::array<double, D> &alpha,
              double coef,
              const Coord<D> &pos = {},
              const std::array<int, D> &power = {});
    GaussPoly(const GaussPoly<D, T> &gp);
    GaussPoly(const GaussFunc<D, T> &gf);
    GaussPoly<D, T> &operator=(const GaussPoly<D, T> &gp) = delete;
    Gaussian<D, T> *copy() const override;
    ~GaussPoly();

    double calcSquareNorm() const override;

    T evalf(const Coord<D> &r) const override;
    T evalf1D(double r, int dim) const override;

    GaussExp<D, T> asGaussExp() const override;
    GaussPoly differentiate(int dir) const override;

    void multInPlace(const GaussPoly<D, T> &rhs);
    void operator*=(const GaussPoly<D, T> &rhs) { multInPlace(rhs); }
    GaussPoly<D, T> mult(const GaussPoly<D, T> &rhs);
    GaussPoly<D, T> mult(double c);
    GaussPoly<D, T> operator*(const GaussPoly<D, T> &rhs) { return mult(rhs); }
    GaussPoly<D, T> operator*(double c) { return mult(c); }

    const Eigen::VectorXd &getPolyCoefs(int i) const { return poly[i]->getCoefs(); }
    Eigen::VectorXd &getPolyCoefs(int i) { return poly[i]->getCoefs(); }
    const Polynomial &getPoly(int i) const { return *poly[i]; }
    Polynomial &getPoly(int i) { return *poly[i]; }

    void setPow(int d, int pow) override;
    void setPow(const std::array<int, D> &pow) override;
    void setPoly(int d, Polynomial &poly);


private:
    Polynomial *poly[D];

    void fillCoefPowVector(std::vector<double> &coefs, std::vector<int *> &power, int pow[D], int dir) const;
    void fillCoefPowVector(std::vector<double> &coefs,
                           std::vector<int *> &power,
                           std::array<int, D> &pow,
                           int dir) const;
    std::ostream &print(std::ostream &o) const override;
};

} // namespace mrcpp
