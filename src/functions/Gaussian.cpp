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
 *
 * \date May 25, 2010
 * \author Stig Rune Jensen
 *		   CTCC, University of Tromsø
 *
 */

#include <numeric>

#include "Gaussian.h"
#include "GaussExp.h"
#include "function_utils.h"
#include "trees/NodeIndex.h"
#include "utils/Printer.h"
#include "utils/details.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

template <int D>
Gaussian<D>::Gaussian(double a, double c, const Coord<D> &r, const std::array<int, D> &p)
        : screen(false)
        , coef(c)
        , power(p)
        , pos(r) {
    this->alpha.fill(a);
}

template <int D>
Gaussian<D>::Gaussian(const std::array<double, D> &a, double c, const Coord<D> &r, const std::array<int, D> &p)
        : screen(false)
        , coef(c)
        , power(p)
        , alpha(a)
        , pos(r) {}

template <int D> void Gaussian<D>::multPureGauss(const Gaussian<D> &lhs, const Gaussian<D> &rhs) {

    auto newAlpha = std::array<double, D>{};
    auto mju = std::array<double, D>{};
    for (auto d = 0; d < D; d++) {
        newAlpha[d] = lhs.alpha[d] + rhs.alpha[d];
        mju[d] = (lhs.alpha[d] * rhs.alpha[d]) / newAlpha[d];
    }
    auto newPos = std::array<double, D>{};
    auto relPos = std::array<double, D>{};

    double newCoef = 1.0;
    for (int d = 0; d < D; d++) {
        newPos[d] = (lhs.alpha[d] * lhs.pos[d] + rhs.alpha[d] * rhs.pos[d]) / newAlpha[d];
        relPos[d] = lhs.pos[d] - rhs.pos[d];
        newCoef *= std::exp(-mju[d] * std::pow(relPos[d], 2.0));
    }
    setExp(newAlpha);
    setPos(newPos);
    setCoef(newCoef);
}

template <int D> void Gaussian<D>::calcScreening(double nStdDev) {
    assert(nStdDev > 0);
    if (not this->isBounded()) {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
    }
    for (int d = 0; d < D; d++) {
        double limit = std::sqrt(nStdDev / this->alpha[d]);
        this->A[d] = this->pos[d] - limit;
        this->B[d] = this->pos[d] + limit;
    }
    screen = true;
}

template <int D> bool Gaussian<D>::checkScreen(int n, const int *l) const {
    if (not getScreen()) { return false; }
    double length = std::pow(2.0, -n);
    const double *A = this->getLowerBounds();
    const double *B = this->getUpperBounds();
    for (int d = 0; d < D; d++) {
        double a = length * l[d];
        double b = length * (l[d] + 1);
        if (a > B[d] or b < A[d]) { return true; }
    }
    return false;
}

template <int D> bool Gaussian<D>::isVisibleAtScale(int scale, int nQuadPts) const {
    for (auto &alp : this->alpha) {
        double stdDeviation = std::pow(2.0 * alp, -0.5);
        auto visibleScale = static_cast<int>(-std::floor(std::log2(nQuadPts * 0.5 * stdDeviation)));

        if (scale < visibleScale) return false;
    }

    return true;
}

template <int D> bool Gaussian<D>::isZeroOnInterval(const double *a, const double *b) const {
    for (int i = 0; i < D; i++) {
        double stdDeviation = std::pow(2.0 * this->alpha[i], -0.5);
        double gaussBoxMin = this->pos[i] - 5.0 * stdDeviation;
        double gaussBoxMax = this->pos[i] + 5.0 * stdDeviation;
        if (a[i] > gaussBoxMax or b[i] < gaussBoxMin) return true;
    }
    return false;
}

template <int D> void Gaussian<D>::evalf(const MatrixXd &points, MatrixXd &values) const {
    assert(points.cols() == D);
    assert(points.cols() == values.cols());
    assert(points.rows() == values.rows());
    for (int d = 0; d < D; d++) {
        for (int i = 0; i < points.rows(); i++) { values(i, d) = evalf1D(points(i, d), d); }
    }
}

template <int D> double Gaussian<D>::getMaximumStandardDiviation() const {

    if (details::are_all_equal<D>(this->getExp())) {
        auto exponent = this->getExp()[0];
        return 1.0 / std::sqrt(2.0 * exponent);
    } else {
        auto exponents = this->getExp();
        auto standard_deviations = std::array<double, D>{};
        for (auto i = 0; i < D; i++) { standard_deviations[i] = 1.0 / std::sqrt(2.0 * exponents[i]); }
        return *std::max_element(standard_deviations.begin(), standard_deviations.end());
    }
}

template <int D> double Gaussian<D>::calcOverlap(const Gaussian<D> &inp) const {
    const auto &bra_exp = this->asGaussExp(); // Make sure all entries are GaussFunc
    const auto &ket_exp = inp.asGaussExp();   // Make sure all entries are GaussFunc

    double S = 0.0;
    for (int i = 0; i < bra_exp.size(); i++) {
        const auto &bra_i = static_cast<const GaussFunc<D> &>(bra_exp.getFunc(i));
        for (int j = 0; j < ket_exp.size(); j++) {
            const auto &ket_j = static_cast<const GaussFunc<D> &>(ket_exp.getFunc(j));
            S += function_utils::calc_overlap(bra_i, ket_j);
        }
    }
    return S;
}

/** @brief Generates a GaussExp that is semi-periodic around a unit-cell
 *
 * @returns Semi-periodic version of a Gaussian around a unit-cell
 * @param[in] period: The period of the unit cell
 * @param[in] nStdDev: Number of standard diviations covered in each direction. Default 4.0
 *
 * @details nStdDev = 1, 2, 3 and 4 ensures atleast 68.27%, 95.45%, 99.73% and 99.99% of the
 * integral is conserved with respect to the integration limits.
 *
 */
template <int D> GaussExp<D> Gaussian<D>::periodify(const std::array<double, D> &period, double nStdDev) const {
    GaussExp<D> gauss_exp;
    auto pos_vec = std::vector<Coord<D>>();

    auto x_std = nStdDev * this->getMaximumStandardDiviation();

    // This lambda function  calculates the number of neighbooring cells
    // requred to keep atleast x_stds of the integral conserved in the
    // unit-cell.
    auto neighbooring_cells = [period, x_std](auto pos) {
        auto needed_cells_vec = std::vector<int>();
        for (auto i = 0; i < D; i++) {
            auto upper_bound = pos[i] + x_std;
            // number of cells upp and down relative to the center of the Gaussian
            needed_cells_vec.push_back(std::ceil(upper_bound / period[i]));
        }

        return *std::max_element(needed_cells_vec.begin(), needed_cells_vec.end());
    };

    // Finding starting position
    auto startpos = this->getPos();

    for (auto d = 0; d < D; d++) {
        startpos[d] = std::fmod(startpos[d], period[d]);
        if (startpos[d] < 0) startpos[d] += period[d];
    }

    auto nr_cells_upp_and_down = neighbooring_cells(startpos);
    for (auto d = 0; d < D; d++) { startpos[d] -= nr_cells_upp_and_down * period[d]; }

    auto tmp_pos = startpos;
    std::vector<double> v(2 * nr_cells_upp_and_down + 1);
    std::iota(v.begin(), v.end(), 0.0);
    auto cart = math_utils::cartesian_product(v, D);
    for (auto &c : cart) {
        for (auto i = 0; i < D; i++) c[i] *= period[i];
    }
    // Shift coordinates
    for (auto &c : cart) std::transform(c.begin(), c.end(), tmp_pos.begin(), c.begin(), std::plus<double>());
    // Go from vector to mrcpp::Coord
    for (auto &c : cart) {
        mrcpp::Coord<D> pos;
        std::copy_n(c.begin(), D, pos.begin());
        pos_vec.push_back(pos);
    }

    for (auto &pos : pos_vec) {
        auto *gauss = this->copy();
        gauss->setPos(pos);
        gauss_exp.append(*gauss);
        delete gauss;
    }

    return gauss_exp;
}

template class Gaussian<1>;
template class Gaussian<2>;
template class Gaussian<3>;

} // namespace mrcpp
