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

    void makePeriodic(const std::array<double, D> &period, double nStdDev = 4.0);

    virtual double evalfCore(const Coord<D> &r) const = 0;

    double evalf(const Coord<D> &r) const;
    virtual double evalf(double r, int dim) const = 0;
    void evalf(const Eigen::MatrixXd &points, Eigen::MatrixXd &values) const;

    virtual double calcSquareNorm() = 0;

    /** @brief Compute analytic overlap between two Gaussians
     *  @param[in] this: Left hand Gaussian
     *  @param[in] rhs: Right hand Gaussian
     *  @returns Overlap integral
     */
    virtual double calcOverlap(GaussFunc<D> &rhs) = 0;

    /** @brief Compute analytic overlap between two Gaussians
     *  @param[in] this: Left hand Gaussian
     *  @param[in] rhs: Right hand Gaussian
     *  @returns Overlap integral
     */
    virtual double calcOverlap(GaussPoly<D> &rhs) = 0;

    /** @brief Compute analytic derivative of Gaussian
     *  @param[in] dir: Cartesian direction of derivative
     *  @returns New GaussPoly
     */
    virtual GaussPoly<D> differentiate(int dir) = 0;

    void calcScreening(double stdDeviations);

    /** @returns Squared L2 norm of function \f$ ||f||^2 \f$ */
    double getSquareNorm() {
        if (this->squareNorm < 0.0) { calcSquareNorm(); }
        return this->squareNorm;
    }
    /** @brief Rescale function by its norm \f$ ||f||^{-1} \f$ */
    void normalize() {
        double norm = std::sqrt(getSquareNorm());
        multConstInPlace(1.0 / norm);
    }
    void multPureGauss(const Gaussian<D> &lhs, const Gaussian<D> &rhs);
    void multConstInPlace(double c) {
        this->coef *= c;
        this->calcSquareNorm();
    }
    void operator*=(double c) { multConstInPlace(c); }

    bool checkScreen(int n, const int *l) const;
    int getPower(int i) const { return power[i]; }
    bool getScreen() const { return screen; }
    const std::array<int, D> &getPower() const { return power; }
    const std::array<double, D> &getPos() const { return pos; }
    double getCoef() const { return coef; }
    std::array<double, D> getExp() const { return alpha; }
    std::array<double, D> getPeriod() const { return period; }

    double getMaximumStandardDiviation() const;

    virtual void setPower(const std::array<int, D> &power) = 0;
    virtual void setPower(int d, int power) = 0;
    void setScreen(bool _screen) { this->screen = _screen; }
    void setCoef(double cf) {
        this->coef = cf;
        this->squareNorm = -1.0;
    }
    void setExp(double _alpha) {
        this->alpha.fill(_alpha);
        this->squareNorm = -1.0;
    }
    void setExp(const std::array<double, D> &_alpha) {
        this->alpha = _alpha;
        this->squareNorm = -1.0;
    }
    void setPos(const std::array<double, D> &r) { this->pos = r; }

    friend std::ostream &operator<<(std::ostream &o, const Gaussian<D> &gauss) { return gauss.print(o); }

    friend class GaussExp<D>;

protected:
    bool screen;
    double coef;                 /**< constant factor */
    std::array<int, D> power;    /**< max power in each dim  */
    std::array<double, D> alpha; /**< exponent  */
    Coord<D> pos;                /**< center  */
    std::array<double, D> period{};
    double squareNorm;

    std::shared_ptr<GaussExp<D>> gauss_exp{nullptr};

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

    virtual std::ostream &print(std::ostream &o) const = 0;
};

} // namespace mrcpp
