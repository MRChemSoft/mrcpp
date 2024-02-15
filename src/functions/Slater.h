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
 *
 *  Base class for Slater type functions
 */

#pragma once

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <memory>

#include "MRCPP/mrcpp_declarations.h"
#include "RepresentableFunction.h"

namespace mrcpp {

template <int D> class Slater : public RepresentableFunction<D> {
public:
    Slater(double a, double c, const Coord<D> &r);
    Slater(const std::array<double, D> &a, double c, const Coord<D> &r);
    Slater<D> &operator=(const Slater<D> &gp) = delete;
    virtual Slater<D> *copy() const = 0;
    virtual ~Slater() = default;

    virtual double evalf(const Coord<D> &r) const = 0;
    virtual double evalf1D(double r, int dim) const = 0;
    void evalf(const Eigen::MatrixXd &points, Eigen::MatrixXd &values) const;

    virtual double calcSquareNorm() const = 0;

    void calcScreening(double stdDeviations);

    /** @brief Rescale function by its norm \f$ ||f||^{-1} \f$ */
    void normalize() {
        double norm = std::sqrt(calcSquareNorm());
        multConstInPlace(1.0 / norm);
    }
    void multConstInPlace(double c) { this->coef *= c; }
    void operator*=(double c) { multConstInPlace(c); }

    bool getScreen() const { return screen; }
    bool checkScreen(int n, const int *l) const;

    double getCoef() const { return coef; }
    double getExp(int i) const { return alpha[i]; }
    const std::array<double, D> &getPos() const { return pos; }
    std::array<double, D> getExp() const { return alpha; }

    void setScreen(bool _screen) { this->screen = _screen; }
    void setCoef(double cf) { this->coef = cf; }
    void setExp(double _alpha) { this->alpha.fill(_alpha); }
    void setExp(const std::array<double, D> &_alpha) { this->alpha = _alpha; }
    void setPos(const std::array<double, D> &r) { this->pos = r; }

    friend std::ostream &operator<<(std::ostream &o, const Slater<D> &gauss) { return gauss.print(o); }

    friend class GaussExp<D>;

protected:
    bool screen;
    double coef;                 /**< constant factor */
    double alpha;                /**< exponent  */
    Coord<D> pos;                /**< center  */

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

    virtual std::ostream &print(std::ostream &o) const = 0;
};

} // namespace mrcpp
