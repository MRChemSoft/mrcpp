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

/** @returns New GaussPoly object
 *  @param[in] beta: Exponent, \f$ e^{-\beta r^2} \f$
 *  @param[in] alpha: Coefficient, \f$ \alpha e^{-r^2} \f$
 *  @param[in] pos: Position \f$ (x - pos[0]), (y - pos[1]), ... \f$
 *  @param[in] pow: Max polynomial degree, \f$ P_0(x), P_1(y), ... \f$
 */
template <int D, typename T>
GaussPoly<D, T>::GaussPoly(double beta, double alpha, const Coord<D> &pos, const std::array<int, D> &power)
        : Gaussian<D, T>(beta, alpha, pos, power) {
    for (auto d = 0; d < D; d++) {
        if (power != std::array<int, D>{}) {
            this->poly[d] = new Polynomial(this->power[d]);
        } else {
            this->poly[d] = nullptr;
        }
    }
}

template <int D, typename T>
GaussPoly<D, T>::GaussPoly(const std::array<double, D> &beta,
                        double alpha,
                        const Coord<D> &pos,
                        const std::array<int, D> &pow)
        : Gaussian<D, T>(beta, alpha, pos, pow) {
    for (auto d = 0; d < D; d++) {
        if (pow != std::array<int, D>{}) {
            this->poly[d] = new Polynomial(this->power[d]);
        } else {
            this->poly[d] = nullptr;
        }
    }
}

template <int D, typename T>
GaussPoly<D, T>::GaussPoly(const GaussPoly<D, T> &gp)
        : Gaussian<D, T>(gp) {
    for (int d = 0; d < D; d++) { poly[d] = new Polynomial(gp.getPoly(d)); }
}

template <int D, typename T>
GaussPoly<D, T>::GaussPoly(const GaussFunc<D, T> &gf)
        : Gaussian<D, T>(gf) {
    for (int d = 0; d < D; d++) {
        int order = this->getPower(d);
        poly[d] = new Polynomial(order);
        VectorXd coefs = VectorXd::Zero(order + 1);
        coefs[order] = 1.0;
        poly[d]->setCoefs(coefs);
        // poly[d]->unsetBounds();
    }
}

template <int D, typename T> GaussPoly<D, T>::~GaussPoly() {
    for (int i = 0; i < D; i++) { delete poly[i]; }
}

template <int D, typename T> Gaussian<D, T> *GaussPoly<D, T>::copy() const {
    auto *gauss = new GaussPoly<D, T>(*this);
    return gauss;
}

template <int D, typename T> double GaussPoly<D, T>::calcSquareNorm() const {
    GaussExp<D, T> this_exp = this->asGaussExp();
    double norm = 0.0;
    for (int i = 0; i < this_exp.size(); i++) {
        auto func_i = static_cast<GaussFunc<D, T> &>(this_exp.getFunc(i));
        for (int j = 0; j < this_exp.size(); j++) {
            auto func_j = static_cast<GaussFunc<D, T> &>(this_exp.getFunc(j));
            norm += function_utils::calc_overlap(func_i, func_j);
        }
    }
    return norm;
}

template <int D, typename T> T GaussPoly<D, T>::evalf(const Coord<D> &r) const {
    if (this->getScreen()) {
        for (int d = 0; d < D; d++) {
            if (r[d] < this->A[d] or r[d] > this->B[d]) { return 0.0; }
        }
    }
    double q2 = 0.0, p2 = 1.0;
    for (int d = 0; d < D; d++) {
        // assert(this->poly[d]->getCheckBounds() == false);
        double q = r[d] - this->pos[d];
        q2 += this->alpha[d] * q * q;
        p2 *= poly[d]->evalf(r[d] - this->pos[d]);
    }
    return this->coef * p2 * std::exp(-q2);
}

template <int D, typename T> T GaussPoly<D, T>::evalf1D(const double r, int d) const {
    // NOTE!
    //     This function evaluation will give the first dimension the full coef
    //     amplitude, leaving all other directions with amplitude 1.0. This is to
    //     avoid expensive d-root evaluation when distributing the amplitude
    //     equally to all dimensions.

    if (this->getScreen()) {
        if ((r < this->A[d]) or (r > this->B[d])) { return 0.0; }
    }
    // assert(this->poly[d]->getCheckBounds() == false);
    double q2 = 0.0, p2 = 1.0;
    double q = (r - this->pos[d]);
    q2 += q * q;
    p2 *= poly[d]->evalf(q);
    if (d == 0) { p2 *= this->coef; }
    return p2 * std::exp(-this->alpha[d] * q2);
}

template <int D, typename T> GaussExp<D, T> GaussPoly<D, T>::asGaussExp() const {
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

    GaussExp<D, T> gexp;
    for (int i = 0; i < nTerms; i++) {
        double coef = coefs[i];
        for (int d = 0; d < D; d++) pow[d] = power[i][d];
        if (coef != 0.0) {
            GaussFunc<D, T> gFunc(alpha, coef, pos, pow);
            gexp.append(gFunc);
        }
    }
    for (auto &i : power) { delete[] i; }
    return gexp;
}

template <int D, typename T> GaussPoly<D, T> GaussPoly<D, T>::differentiate(int dir) const {
    NOT_IMPLEMENTED_ABORT;
}

template <int D, typename T> void GaussPoly<D, T>::multInPlace(const GaussPoly<D, T> &rhs) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D, typename T>
void GaussPoly<D, T>::fillCoefPowVector(std::vector<double> &coefs, std::vector<int *> &power, int pow[D], int dir) const {
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

template <int D, typename T>
void GaussPoly<D, T>::fillCoefPowVector(std::vector<double> &coefs,
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

template <int D, typename T> GaussPoly<D, T> GaussPoly<D, T>::mult(const GaussPoly<D, T> &rhs) {
    NOT_IMPLEMENTED_ABORT;
    /*
    GaussPoly<D, T> &lhs = *this;
    GaussPoly<D, T> result;
    result.multPureGauss(lhs, rhs);
    for (int d = 0; d < D; d++) {
        double newPos = result.getPos()[d];
        int lhsPow = lhs.getPower(d);
        Polynomial lhsPoly(lhsPow);
        lhsPoly.clearCoefs();
        for (int p = 0; p <= lhsPow; p++) {
            Polynomial tmpPoly(newPos - lhs.getPos()[d], p);
            tmpPoly *= lhs.getPolyCoefs(d)[p];
            lhsPoly += tmpPoly;
        }

        int rhsPow = rhs.getPower(d);
        Polynomial rhsPoly(rhsPow);
        rhsPoly.clearCoefs();
        for (int p = 0; p <= rhsPow; p++) {
            Polynomial tmpPoly(newPos - rhs.getPos()[d], p);
            tmpPoly *= rhs.getPolyCoefs(d)[p];
            rhsPoly += tmpPoly;
        }
        Polynomial newPoly = lhsPoly * rhsPoly;
        result.setPoly(d, newPoly);
    }
    result.setCoef(result.getCoef() * lhs.getCoef() * rhs.getCoef());
    return result;
    */
}

/** @brief Multiply GaussPoly by scalar
 *  @param[in] c: Scalar to multiply
 *  @returns New GaussPoly
 */
template <int D, typename T> GaussPoly<D, T> GaussPoly<D, T>::mult(double c) {
    GaussPoly<D, T> g = *this;
    g.coef *= c;
    return g;
}

template <int D, typename T> void GaussPoly<D, T>::setPow(int d, int pow) {
    if (poly[d] != nullptr) { delete poly[d]; }
    poly[d] = new Polynomial(pow);
}

template <int D, typename T> void GaussPoly<D, T>::setPow(const std::array<int, D> &pow) {
    for (int d = 0; d < D; d++) {
        if (poly[d] != nullptr) { delete poly[d]; }
        poly[d] = new Polynomial(pow[d]);
    }
}

/** @brief Set polynomial in given dimension
 *
 *  @param[in] d: Cartesian direction
 *  @param[in] poly: Polynomial to set
 */
template <int D, typename T> void GaussPoly<D, T>::setPoly(int d, Polynomial &poly) {
    if (this->poly[d] != nullptr) { delete this->poly[d]; }
    this->poly[d] = new Polynomial(poly);
    this->power[d] = poly.getOrder();
}

template <int D, typename T> std::ostream &GaussPoly<D, T>::print(std::ostream &o) const {
    auto is_array = details::are_all_equal<D>(this->getExp());

    // If all of the values in the exponential are the same only
    // one is printed, else, all of them are printed
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

template class GaussPoly<1, double>;
template class GaussPoly<2, double>;
template class GaussPoly<3, double>;

template class GaussPoly<1, ComplexDouble>;
template class GaussPoly<2, ComplexDouble>;
template class GaussPoly<3, ComplexDouble>;

} // namespace mrcpp
