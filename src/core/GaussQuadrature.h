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

#include "functions/RepresentableFunction.h"

namespace mrcpp {

const int MaxGaussOrder = 42;
static const double EPS = 3.0e-12;
static const int NewtonMaxIter = 10;
static const int MaxQuadratureDim = 7;

class GaussQuadrature final {
public:
    GaussQuadrature(int k, double a = -1.0, double b = 1.0, int inter = 1);

    double integrate(RepresentableFunction<1> &func) const;
    double integrate(RepresentableFunction<2> &func) const;
    double integrate(RepresentableFunction<3> &func) const;

    void setIntervals(int i);
    void setBounds(double a, double b);

    int getIntervals() const { return this->intervals; }
    double getUpperBound() const { return this->B; }
    double getLowerBound() const { return this->A; }

    const Eigen::VectorXd &getRoots() const { return this->roots; }
    const Eigen::VectorXd &getWeights() const { return this->weights; }
    const Eigen::VectorXd &getUnscaledRoots() const { return this->unscaledRoots; }
    const Eigen::VectorXd &getUnscaledWeights() const { return this->unscaledWeights; }

protected:
    int order;
    double A;
    double B;
    int intervals;
    int npts;
    Eigen::VectorXd roots;
    Eigen::VectorXd weights;
    Eigen::VectorXd unscaledRoots;
    Eigen::VectorXd unscaledWeights;

    void rescaleRoots(Eigen::VectorXd &rts, double a, double b, int inter = 1) const;
    void rescaleWeights(Eigen::VectorXd &wgts, double a, double b, int inter = 1) const;

    void calcScaledPtsWgts();
    int calcGaussPtsWgts();

    double integrate_nd(RepresentableFunction<3> &func, int axis = 0) const;
};

} // namespace mrcpp
