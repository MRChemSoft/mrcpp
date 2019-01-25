/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
 *
 *
 * \date May 25, 2010
 * \author Stig Rune Jensen
 *		   CTCC, University of Tromsø
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include <vector>

#include "Gaussian.h"
#include "Polynomial.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

template <int D> class GaussPoly final : public Gaussian<D> {
public:
    GaussPoly(double alpha = 0.0, double coef = 1.0, const Coord<D> &pos = {}, const std::array<int, D> &power = {});
    GaussPoly(const std::array<double, D> &alpha,
              double coef,
              const Coord<D> &pos = {},
              const std::array<int, D> &power = {});
    GaussPoly(const GaussPoly<D> &gp);
    GaussPoly(const GaussFunc<D> &gf);
    GaussPoly<D> &operator=(const GaussPoly<D> &gp) = delete;
    Gaussian<D> *copy() const;
    ~GaussPoly();

    double calcSquareNorm();

    double evalf(const Coord<D> &r) const;
    double evalf(double r, int dim) const;

    double calcOverlap(GaussFunc<D> &b);
    double calcOverlap(GaussPoly<D> &b);

    GaussPoly differentiate(int dir);

    void multInPlace(const GaussPoly<D> &rhs);
    void operator*=(const GaussPoly<D> &rhs) { multInPlace(rhs); }
    GaussPoly<D> mult(const GaussPoly<D> &rhs);
    GaussPoly<D> mult(double c);
    GaussPoly<D> operator*(const GaussPoly<D> &rhs) { return mult(rhs); }
    GaussPoly<D> operator*(double c) { return mult(c); }

    const Eigen::VectorXd &getPolyCoefs(int i) const { return poly[i]->getCoefs(); }
    Eigen::VectorXd &getPolyCoefs(int i) { return poly[i]->getCoefs(); }
    const Polynomial &getPoly(int i) const { return *poly[i]; }
    Polynomial &getPoly(int i) { return *poly[i]; }

    void setPower(int d, int pow);
    void setPower(const std::array<int, D> &pow);
    // void setPower(const int pow[D]);
    void setPoly(int d, Polynomial &poly);

    void fillCoefPowVector(std::vector<double> &coefs, std::vector<int *> &power, int pow[D], int dir) const;
    void fillCoefPowVector(std::vector<double> &coefs,
                           std::vector<int *> &power,
                           std::array<int, D> &pow,
                           int dir) const;

private:
    Polynomial *poly[D];

    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp
