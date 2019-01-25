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
 * \breif Tools to handle gaussian functions and/or expansions in gaussian functions.
 * Implemented as a template class in the dimensionality:
 *
 *- Monodimensional gaussian (Gaussian<1>):
 *
 * \f$ g(x) = c e^{-\alpha (x-x^0)^2} \f$
 *
 * - Multidimensional gaussian (Gaussian<d>):
 *
 * \f$ G(x) = \prod_{i=1}^d g_i(x_i)
 * = \prod_{i=1}^d \c_i e^{-\alpha_i (x_i-x_i^0)^2} \f$
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "Gaussian.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

template <int D> class GaussFunc final : public Gaussian<D> {
public:
    GaussFunc(double alpha, double coef, const Coord<D> &pos = {}, const std::array<int, D> &pow = {})
            : Gaussian<D>(alpha, coef, pos, pow) {}
    GaussFunc(const std::array<double, D> &alpha, double coef, const Coord<D> &pos, const std::array<int, D> &pow)
            : Gaussian<D>(alpha, coef, pos, pow) {}
    GaussFunc(const GaussFunc<D> &gf)
            : Gaussian<D>(gf) {}
    GaussFunc<D> &operator=(const GaussFunc<D> &gp) = delete;
    Gaussian<D> *copy() const;

    double calcCoulombEnergy(GaussFunc<D> &gf);
    double calcSquareNorm();

    double evalf(const Coord<D> &r) const;
    double evalf(double r, int dim) const;

    static double calcOverlap(GaussFunc<D> &a, GaussFunc<D> &b);
    double calcOverlap(GaussFunc<D> &b);
    double calcOverlap(GaussPoly<D> &b);

    GaussPoly<D> differentiate(int dir);

    void multInPlace(const GaussFunc<D> &g);
    void operator*=(const GaussFunc<D> &gf) { multInPlace(gf); }
    GaussPoly<D> mult(const GaussFunc<D> &g);
    GaussFunc<D> mult(double d);
    GaussPoly<D> operator*(const GaussFunc<D> &g) { return this->mult(g); }
    GaussFunc<D> operator*(double d) { return this->mult(d); }

    void setPower(int d, int power) {
        this->power[d] = power;
        this->squareNorm = -1.0;
    }
    void setPower(const std::array<int, D> &power) {
        this->power = power;
        this->squareNorm = -1.0;
    }

private:
    static double ObaraSaika_ab(int power_a, int power_b, double pos_a, double pos_b, double expo_a, double expo_b);

    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp
