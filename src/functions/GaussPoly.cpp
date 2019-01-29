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
 *
 * \date May 26, 2010
 * \author Stig Rune Jensen
 *		   CTCC, University of Tromsø
 *
 */

#include "vector"

#include "GaussExp.h"
#include "GaussFunc.h"
#include "GaussPoly.h"
#include "utils/Printer.h"
#include "utils/details.h"

using namespace Eigen;

namespace mrcpp {

template <int D>
GaussPoly<D>::GaussPoly(double alpha, double coef, const Coord<D> &pos, const std::array<int, D> &power)
        : Gaussian<D>(alpha, coef, pos, power) {
    for (auto d = 0; d < D; d++) {
        if (power != std::array<int, D>{}) {
            this->poly[d] = new Polynomial(this->power[d]);
            // this->poly[d]->unsetBounds();
        } else {
            this->poly[d] = nullptr;
        }
    }
}

template <int D>
GaussPoly<D>::GaussPoly(const std::array<double, D> &alpha,
                        double coef,
                        const Coord<D> &pos,
                        const std::array<int, D> &power)
        : Gaussian<D>(alpha, coef, pos, power) {
    for (auto d = 0; d < D; d++) {
        if (power != std::array<int, D>{}) {
            this->poly[d] = new Polynomial(this->power[d]);
            // this->poly[d]->unsetBounds();
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
        // poly[d]->unsetBounds();
    }
}

template <int D> GaussPoly<D>::~GaussPoly() {
    for (int i = 0; i < D; i++) { delete poly[i]; }
}

template <int D> Gaussian<D> *GaussPoly<D>::copy() const {
    auto *gauss = new GaussPoly<D>(*this);
    return gauss;
}

template <int D> double GaussPoly<D>::calcOverlap(GaussPoly<D> &b) {
    GaussExp<D> gExp(*this);
    GaussExp<D> fExp(b);

    double overlap = 0.0;
    for (int i = 0; i < fExp.size(); i++) {
        auto &fFunc = static_cast<GaussFunc<D> &>(fExp.getFunc(i));
        for (int j = 0; j < gExp.size(); j++) { overlap += gExp.getFunc(j).calcOverlap(fFunc); }
    }
    return overlap;
}

template <int D> double GaussPoly<D>::calcOverlap(GaussFunc<D> &b) {
    return b.calcOverlap(*this);
}

template <int D> double GaussPoly<D>::calcSquareNorm() {
    this->squareNorm = this->calcOverlap(*this);
    return this->squareNorm;
}

template <int D> double GaussPoly<D>::evalf(const Coord<D> &r) const {
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

/** NOTE!
 *	This function evaluation will give the first dimension the full coef
 *	amplitude, leaving all other directions with amplitude 1.0. This is to
 *	avoid expensive d-root evaluation when distributing the amplitude
 *	equally to all dimensions.
 */
template <int D> double GaussPoly<D>::evalf(const double r, int d) const {
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

template <int D> GaussPoly<D> GaussPoly<D>::differentiate(int dir) {
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
            int *newPow = new int[D];
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
            int *newPow = new int[D];
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
    /*
    GaussPoly<D> &lhs = *this;
    GaussPoly<D> result;
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

template <int D> GaussPoly<D> GaussPoly<D>::mult(double c) {
    GaussPoly<D> g = *this;
    g.coef *= c;
    return g;
}

template <int D> void GaussPoly<D>::setPower(int d, int pow) {
    if (poly[d] != nullptr) { delete poly[d]; }
    poly[d] = new Polynomial(pow);
    this->squareNorm = -1.0;
}

template <int D> void GaussPoly<D>::setPower(const std::array<int, D> &pow) {
    for (int d = 0; d < D; d++) {
        if (poly[d] != nullptr) { delete poly[d]; }
        poly[d] = new Polynomial(pow[d]);
    }
    this->squareNorm = -1.0;
}

template <int D> void GaussPoly<D>::setPoly(int d, Polynomial &poly) {
    if (this->poly[d] != nullptr) { delete this->poly[d]; }
    this->poly[d] = new Polynomial(poly);
    // this->poly[d]->unsetBounds();
    this->power[d] = poly.getOrder();
}

template <int D> std::ostream &GaussPoly<D>::print(std::ostream &o) const {
    auto expTmp = this->getExp();
    auto is_array = details::are_all_equal<D>(expTmp);

    // If all of the values in the exponential are the same only
    // one is printed, else, all of them are printed
    if (!is_array) {
        o << "Exp:   ";
        for (auto &alpha : expTmp) { o << alpha << " "; }
    } else {
        o << "Exp:   " << expTmp[0] << std::endl;
    }
    o << "Coef: " << this->getCoef() << std::endl;
    o << "Pos:   ";
    for (int i = 0; i < D; i++) { o << this->getPos()[i] << " "; }
    o << std::endl;
    for (int i = 0; i < D; i++) {
        o << "Dim " << i << ": order " << this->getPower(i) << std::endl;
        o << this->getPolyCoefs(i) << std::endl;
    }
    return o;
}

template class GaussPoly<1>;
template class GaussPoly<2>;
template class GaussPoly<3>;

} // namespace mrcpp
