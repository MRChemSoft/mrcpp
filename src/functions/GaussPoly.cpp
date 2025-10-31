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

#include "GaussPoly.h"

#include <vector>

#include "function_utils.h"
#include "GaussExp.h"
#include "GaussFunc.h"
#include "utils/Printer.h"
#include "utils/details.h"

using namespace Eigen;

namespace mrcpp {

template <int D>
GaussPoly<D>::GaussPoly(double beta, double alpha, const Coord<D> &pos, const std::array<int, D> &power)
        : Gaussian<D>(beta, alpha, pos, power) {
    for (auto d = 0; d < D; d++) {
        if (power != std::array<int, D>{}) {
            this->poly[d] = new Polynomial(this->power[d]);
        } else {
            this->poly[d] = nullptr;
        }
    }
}

template <int D>
GaussPoly<D>::GaussPoly(const std::array<double, D> &beta,
                        double alpha,
                        const Coord<D> &pos,
                        const std::array<int, D> &pow)
        : Gaussian<D>(beta, alpha, pos, pow) {
    for (auto d = 0; d < D; d++) {
        if (pow != std::array<int, D>{}) {
            this->poly[d] = new Polynomial(this->power[d]);
        } else {
            this->poly[d] = nullptr;
        }
    }
}

template <int D>
GaussPoly<D>::GaussPoly(const GaussPoly<D> &gp)
        : Gaussian<D>(gp) {
    for (int d = 0; d < D; d++) { poly[d] = new Polynomial(gp.getPoly(d)); }
}

template <int D>
GaussPoly<D>::GaussPoly(const GaussFunc<D> &gf)
        : Gaussian<D>(gf) {
    for (int d = 0; d < D; d++) {
        int order = this->getPower(d);
        poly[d] = new Polynomial(order);
        VectorXd coefs = VectorXd::Zero(order + 1);
        coefs[order] = 1.0;
        poly[d]->setCoefs(coefs);
    }
}

template <int D> GaussPoly<D>::~GaussPoly() {
    for (int i = 0; i < D; i++) { delete poly[i]; }
}

template <int D> Gaussian<D> *GaussPoly<D>::copy() const {
    auto *gauss = new GaussPoly<D>(*this);
    return gauss;
}

template<int D> double GaussPoly<D>::calcSquareNorm() const {
    GaussExp<D> this_exp = this->asGaussExp();
    double norm = 0.0;
    for (int i = 0; i < this_exp.size(); i++) {
        auto func_i = static_cast<GaussFunc<D> &>(this_exp.getFunc(i));
        for (int j = 0; j < this_exp.size(); j++) {
            auto func_j = static_cast<GaussFunc<D> &>(this_exp.getFunc(j));
            norm += function_utils::calc_overlap(func_i, func_j);
        }
    }
    return norm;
}

template <int D> double GaussPoly<D>::evalf(const Coord<D> &r) const {
    if (this->getScreen()) {
        for (int d = 0; d < D; d++) {
            if (r[d] < this->A[d] or r[d] > this->B[d]) { return 0.0; }
        }
    }
    double q2 = 0.0, p2 = 1.0;
    for (int d = 0; d < D; d++) {
        double q = r[d] - this->pos[d];
        q2 += this->alpha[d] * q * q;
        p2 *= poly[d]->evalf(r[d] - this->pos[d]);
    }
    return this->coef * p2 * std::exp(-q2);
}

template <int D> double GaussPoly<D>::evalf1D(const double r, int d) const {
    if (this->getScreen()) {
        if ((r < this->A[d]) or (r > this->B[d])) { return 0.0; }
    }
    double q2 = 0.0, p2 = 1.0;
    double q = (r - this->pos[d]);
    q2 += q * q;
    p2 *= poly[d]->evalf(q);
    if (d == 0) { p2 *= this->coef; }
    return p2 * std::exp(-this->alpha[d] * q2);
}

template <int D> GaussExp<D> GaussPoly<D>::asGaussExp() const {
    std::array<int, D> pow;
    std::array<double, D> pos;
    auto alpha = this->getExp();

    int nTerms = 1;
    for (int d = 0; d < D; d++) {
        nTerms *= (this->getPower(d) + 1);
        pos[d] = this->getPos()[d];
    }

    std::vector<double> coefs;
    std::vector<int *> power;

    fillCoefPowVector(coefs, power, pow, D);

    GaussExp<D> gexp;
    for (int i = 0; i < nTerms; i++) {
        double coef = coefs[i];
        for (int d = 0; d < D; d++) pow[d] = power[i][d];
        if (coef != 0.0) {
            GaussFunc<D> gFunc(alpha, coef, pos, pow);
            gexp.append(gFunc);
        }
    }
    for (auto &i : power) { delete[] i; }
    return gexp;
}

template <int D> GaussPoly<D> GaussPoly<D>::differentiate(int dir) const {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> void GaussPoly<D>::multInPlace(const GaussPoly<D> &rhs) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D>
void GaussPoly<D>::fillCoefPowVector(std::vector<double> &coefs, std::vector<int *> &power, int pow[D], int dir) const {
    dir--;
    for (int i = 0; i < this->getPower(dir) + 1; i++) {
        pow[dir] = i;
        if (dir > 0) {
            fillCoefPowVector(coefs, power, pow, dir);
        } else {
            auto *newPow = new int[D];
            double coef = 1.0;
            for (int d = 0; d < D; d++) {
                newPow[d] = pow[d];
                coef *= this->getPolyCoefs(d)[pow[d]];
            }
            coef *= this->getCoef();
            power.push_back(newPow);
            coefs.push_back(coef);
        }
    }
}

template <int D>
void GaussPoly<D>::fillCoefPowVector(std::vector<double> &coefs,
                                     std::vector<int *> &power,
                                     std::array<int, D> &pow,
                                     int dir) const {
    dir--;
    for (int i = 0; i < this->getPower(dir) + 1; i++) {
        pow[dir] = i;
        if (dir > 0) {
            fillCoefPowVector(coefs, power, pow, dir);
        } else {
            auto *newPow = new int[D];
            double coef = 1.0;
            for (int d = 0; d < D; d++) {
                newPow[d] = pow[d];
                coef *= this->getPolyCoefs(d)[pow[d]];
            }
            coef *= this->getCoef();
            power.push_back(newPow);
            coefs.push_back(coef);
        }
    }
}

template <int D> GaussPoly<D> GaussPoly<D>::mult(const GaussPoly<D> &rhs) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> GaussPoly<D> GaussPoly<D>::mult(double c) {
    GaussPoly<D> g = *this;
    g.coef *= c;
    return g;
}

template <int D> void GaussPoly<D>::setPow(int d, int pow) {
    if (poly[d] != nullptr) { delete poly[d]; }
    poly[d] = new Polynomial(pow);
}

template <int D> void GaussPoly<D>::setPow(const std::array<int, D> &pow) {
    for (int d = 0; d < D; d++) {
        if (poly[d] != nullptr) { delete poly[d]; }
        poly[d] = new Polynomial(pow[d]);
    }
}

template <int D> void GaussPoly<D>::setPoly(int d, Polynomial &poly) {
    if (this->poly[d] != nullptr) { delete this->poly[d]; }
    this->poly[d] = new Polynomial(poly);
    this->power[d] = poly.getOrder();
}

template <int D> std::ostream &GaussPoly<D>::print(std::ostream &o) const {
    auto is_array = details::are_all_equal<D>(this->getExp());
    o << "Coef    : " << this->getCoef() << std::endl;
    if (!is_array) {
        o << "Exp     : ";
        for (auto &alpha : this->getExp()) o << alpha << " ";
    } else {
        o << "Exp     : " << this->getExp(0) << std::endl;
    }
    o << "Pos     : ";
    for (int i = 0; i < D; i++) o << this->getPos()[i] << " ";
    o << std::endl;
    o << "Pow     : ";
    for (int i = 0; i < D; i++) o << this->getPower()[i] << " ";
    o << std::endl;
    for (int i = 0; i < D; i++) o << "Poly[" << i << "] : " << this->getPolyCoefs(i).transpose() << std::endl;
    return o;
}

template class GaussPoly<1>;
template class GaussPoly<2>;
template class GaussPoly<3>;

} // namespace mrcpp