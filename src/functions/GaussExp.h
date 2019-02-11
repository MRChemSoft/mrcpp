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

/**
 *
 * - Monodimensional gaussian expansion:
 *
 * \f$ g(x) = \sum_{i=1}^M g_i(x_i)
 * = \sum_{i=1}^M \alpha_i e^{-c_i (x-x^0)^2} \f$
 *
 * - Multidimensional gaussian expansion:
 *
 * \f$ G(x) = \sum_{j=1}^M G_j(x)
 * = \sum_{j=1}^M \prod_{i=1}^d g_ij(x_i)
 * = \sum_{j=1}^M \prod_{i=1}^d \alpha_{ij} e^{-c_{ij} (x_i-x_i^0)^2} \f$
 *
 *
 */

#pragma once

#include <iostream>
#include <vector>

#include "GaussFunc.h"
#include "RepresentableFunction.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

#define GAUSS_EXP_PREC 1.e-10

template <int D> class GaussExp : public RepresentableFunction<D> {
public:
    GaussExp(int nTerms = 0, double prec = GAUSS_EXP_PREC);
    GaussExp(const GaussExp<D> &gExp);
    GaussExp(const GaussPoly<D> &gPoly);
    GaussExp &operator=(const GaussExp<D> &gExp);
    virtual ~GaussExp();

    double calcCoulombEnergy();
    double calcSquareNorm();
    void normalize();

    void calcScreening(double nStdDev = defaultScreening);
    bool isVisibleAtScale(int scale, int nPts) const;
    bool isZeroOnInterval(const double *lb, const double *ub) const;

    double evalf(const Coord<D> &r) const;

    GaussExp<D> differentiate(int dir);

    GaussExp<D> add(GaussExp<D> &g);
    GaussExp<D> add(Gaussian<D> &g);
    GaussExp<D> mult(GaussExp<D> &g);
    GaussExp<D> mult(GaussFunc<D> &g);
    GaussExp<D> mult(GaussPoly<D> &g);
    GaussExp<D> mult(double d);
    void multInPlace(double d);

    GaussExp<D> operator+(GaussExp<D> &g) { return this->add(g); }
    GaussExp<D> operator+(Gaussian<D> &g) { return this->add(g); }
    GaussExp<D> operator*(GaussExp<D> &g) { return this->mult(g); }
    GaussExp<D> operator*(GaussFunc<D> &g) { return this->mult(g); }
    GaussExp<D> operator*(GaussPoly<D> &g) { return this->mult(g); }
    GaussExp<D> operator*(double d) { return this->mult(d); }
    void operator*=(double d) { this->multInPlace(d); }

    double getScreening() const { return screening; }
    std::array<double, D> getExp(int i) const { return this->funcs[i]->getExp(); }
    double getCoef(int i) const { return this->funcs[i]->getCoef(); }
    const std::array<int, D> &getPower(int i) const { return this->funcs[i]->getPower(); }
    const std::array<double, D> &getPos(int i) const { return this->funcs[i]->getPos(); }

    double getSquareNorm() {
        if (squareNorm < 0) { calcSquareNorm(); }
        return squareNorm;
    }

    int size() const { return this->funcs.size(); }
    Gaussian<D> &getFunc(int i) { return *this->funcs[i]; }
    const Gaussian<D> &getFunc(int i) const { return *this->funcs[i]; }

    Gaussian<D> *operator[](int i) { return this->funcs[i]; }
    const Gaussian<D> *operator[](int i) const { return this->funcs[i]; }

    void setFunc(int i, const GaussPoly<D> &g, double c = 1.0);
    void setFunc(int i, const GaussFunc<D> &g, double c = 1.0);

    void setDefaultScreening(double screen);
    void setScreen(bool screen);
    void setExp(int i, double a) { this->funcs[i]->setExp(a); }
    void setCoef(int i, double b) { this->funcs[i]->setCoef(b); }
    void setPower(int i, const std::array<int, D> &power) { this->funcs[i]->setPower(power); }
    void setPos(int i, const std::array<double, D> &pos) { this->funcs[i]->setPos(pos); }

    void append(const Gaussian<D> &g);
    void append(const GaussExp<D> &g);

    friend std::ostream &operator<<(std::ostream &o, const GaussExp<D> &gExp) { return gExp.print(o); }

protected:
    std::vector<Gaussian<D> *> funcs;
    static double defaultScreening;
    double screening{0.0};
    double squareNorm{-1.0};

    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp
