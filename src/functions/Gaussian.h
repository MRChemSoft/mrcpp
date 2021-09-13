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
 *  Base class for Gaussian type functions
 */

#pragma once

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <memory>

#include "MRCPP/mrcpp_declarations.h"
#include "RepresentableFunction.h"

namespace mrcpp {

template <int D> class Gaussian : public RepresentableFunction<D> {
public:
    Gaussian(double a, double c, const Coord<D> &r, const std::array<int, D> &p);
    Gaussian(const std::array<double, D> &a, double c, const Coord<D> &r, const std::array<int, D> &p);
    Gaussian<D> &operator=(const Gaussian<D> &gp) = delete;
    virtual Gaussian<D> *copy() const = 0;
    virtual ~Gaussian() = default;

    virtual double evalf(const Coord<D> &r) const = 0;
    virtual double evalf1D(double r, int dim) const = 0;
    void evalf(const Eigen::MatrixXd &points, Eigen::MatrixXd &values) const;

    double calcOverlap(const Gaussian<D> &inp) const;
    virtual double calcSquareNorm() const = 0;
    virtual GaussExp<D> asGaussExp() const = 0;
    GaussExp<D> periodify(const std::array<double, D> &period, double nStdDev = 4.0) const;

    /** @brief Compute analytic derivative of Gaussian
     *  @param[in] dir: Cartesian direction of derivative
     *  @returns New GaussPoly
     */
    virtual GaussPoly<D> differentiate(int dir) const = 0;

    void calcScreening(double stdDeviations);

    /** @brief Rescale function by its norm \f$ ||f||^{-1} \f$ */
    void normalize() {
        double norm = std::sqrt(calcSquareNorm());
        multConstInPlace(1.0 / norm);
    }
    void multPureGauss(const Gaussian<D> &lhs, const Gaussian<D> &rhs);
    void multConstInPlace(double c) { this->coef *= c; }
    void operator*=(double c) { multConstInPlace(c); }

    bool getScreen() const { return screen; }
    bool checkScreen(int n, const int *l) const;

    int getPower(int i) const { return power[i]; }
    double getCoef() const { return coef; }
    double getExp(int i) const { return alpha[i]; }
    const std::array<int, D> &getPower() const { return power; }
    const std::array<double, D> &getPos() const { return pos; }
    std::array<double, D> getExp() const { return alpha; }

    virtual void setPow(const std::array<int, D> &power) = 0;
    virtual void setPow(int d, int power) = 0;
    void setScreen(bool _screen) { this->screen = _screen; }
    void setCoef(double cf) { this->coef = cf; }
    void setExp(double _alpha) { this->alpha.fill(_alpha); }
    void setExp(const std::array<double, D> &_alpha) { this->alpha = _alpha; }
    void setPos(const std::array<double, D> &r) { this->pos = r; }

    friend std::ostream &operator<<(std::ostream &o, const Gaussian<D> &gauss) { return gauss.print(o); }

    friend class GaussExp<D>;

protected:
    bool screen;
    double coef;                 /**< constant factor */
    std::array<int, D> power;    /**< max power in each dim  */
    std::array<double, D> alpha; /**< exponent  */
    Coord<D> pos;                /**< center  */

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

    double getMaximumStandardDiviation() const;

    virtual std::ostream &print(std::ostream &o) const = 0;
};

} // namespace mrcpp
