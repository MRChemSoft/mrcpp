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

#include "GaussExp.h"

#include <cstdlib>
#include <iostream>

#include "function_utils.h"
#include "GaussFunc.h"
#include "GaussPoly.h"
#include "Gaussian.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

template <int D> double GaussExp<D>::defaultScreening = 10.0;

template <int D> GaussExp<D>::GaussExp(int nTerms) {
    for (int i = 0; i < nTerms; i++) { this->funcs.push_back(nullptr); }
}

template <int D> GaussExp<D>::GaussExp(const GaussExp<D> &gexp) {
    screening = gexp.screening;
    for (unsigned int i = 0; i < gexp.size(); i++) {
        Gaussian<D> *gauss = gexp.funcs[i]->copy();
        this->funcs.push_back(gauss);
    }
}

template <int D> GaussExp<D>::~GaussExp() {
    for (int i = 0; i < size(); i++) {
        if (this->funcs[i] != nullptr) {
            delete this->funcs[i];
            this->funcs[i] = nullptr;
        }
    }
}

template <int D> GaussExp<D> &GaussExp<D>::operator=(const GaussExp<D> &gexp) {
    if (&gexp == this) return *this;
    // screening = gexp.screening;
    this->funcs.clear();
    for (unsigned int i = 0; i < gexp.size(); i++) {
        if (gexp.funcs[i] == nullptr) {
            this->funcs.push_back(nullptr);
        } else {
            Gaussian<D> *gauss = gexp.getFunc(i).copy();
            this->funcs.push_back(gauss);
        }
    }
    return *this;
}

template <int D> double GaussExp<D>::evalf(const Coord<D> &r) const {
    double val = 0.0;
    for (int i = 0; i < this->size(); i++) { val += this->getFunc(i).evalf(r); }
    return val;
}

template <int D> bool GaussExp<D>::isVisibleAtScale(int scale, int nPts) const {
    for (unsigned int i = 0; i < this->size(); i++) {
        if (not this->getFunc(i).isVisibleAtScale(scale, nPts)) { return false; }
    }
    return true;
}

template <int D> bool GaussExp<D>::isZeroOnInterval(const double *lb, const double *ub) const {
    for (unsigned int i = 0; i < this->size(); i++) {
        if (not this->getFunc(i).isZeroOnInterval(lb, ub)) { return false; }
    }
    return true;
}

template <int D> void GaussExp<D>::setFunc(int i, const GaussPoly<D> &g, double c) {
    if (i < 0 or i > (this->size() - 1)) {
        MSG_ERROR("Index out of bounds!");
        return;
    }
    if (this->funcs[i] != nullptr) { delete this->funcs[i]; }
    this->funcs[i] = new GaussPoly<D>(g);
    double coef = this->funcs[i]->getCoef();
    this->funcs[i]->setCoef(c * coef);
}

template <int D> void GaussExp<D>::setFunc(int i, const GaussFunc<D> &g, double c) {
    if (i < 0 or i > (this->size() - 1)) {
        MSG_ERROR("Index out of bounds!");
        return;
    }
    if (this->funcs[i] != nullptr) { delete this->funcs[i]; }
    this->funcs[i] = new GaussFunc<D>(g);
    double coef = this->funcs[i]->getCoef();
    this->funcs[i]->setCoef(c * coef);
}

template <int D> void GaussExp<D>::append(const Gaussian<D> &g) {
    Gaussian<D> *gp = g.copy();
    this->funcs.push_back(gp);
}

template <int D> void GaussExp<D>::append(const GaussExp<D> &g) {
    for (int i = 0; i < g.size(); i++) {
        Gaussian<D> *gp = g.getFunc(i).copy();
        this->funcs.push_back(gp);
    }
}

template <int D> GaussExp<D> GaussExp<D>::differentiate(int dir) const {
    assert(dir >= 0 and dir < D);
    GaussExp<D> result;
    for (int i = 0; i < this->size(); i++) result.append(this->getFunc(i).differentiate(dir));
    return result;
}

template <int D> GaussExp<D> GaussExp<D>::add(GaussExp<D> &g) {
    int nsum = this->size() + g.size();
    GaussExp<D> sum = GaussExp<D>(nsum);

    int n = 0;
    for (int i = 0; i < this->size(); i++) {
        sum.funcs[n] = this->funcs[i]->copy();
        n++;
    }
    for (int i = 0; i < g.size(); i++) {
        sum.funcs[n] = g.funcs[i]->copy();
        n++;
    }

    return sum;
}

template <int D> GaussExp<D> GaussExp<D>::add(Gaussian<D> &g) {
    int nsum = this->size() + 1;
    GaussExp<D> sum = GaussExp<D>(nsum);
    for (int n = 0; n < this->size(); n++) { sum.funcs[n] = this->getFunc(n).copy(); }
    sum.funcs[this->size()] = g.copy();
    return sum;
}

template <int D> GaussExp<D> GaussExp<D>::mult(GaussExp<D> &gexp) {
    GaussExp<D> result;
    for (int i = 0; i < this->size(); i++) {
        for (int j = 0; j < gexp.size(); j++) {
            if (auto *f = dynamic_cast<GaussFunc<D> *>(this->funcs[i])) {
                if (auto *g = dynamic_cast<GaussFunc<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*g) * (*f);
                    result.append(newTerm);
                } else if (auto *g = dynamic_cast<GaussPoly<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*g) * (*f);
                    result.append(newTerm);
                } else {
                    MSG_ABORT("Invalid Gaussian type!");
                }
            } else if (auto *f = dynamic_cast<GaussPoly<D> *>(this->funcs[i])) {
                if (auto *g = dynamic_cast<GaussFunc<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*f) * (*g);
                    result.append(newTerm);
                } else if (auto *g = dynamic_cast<GaussPoly<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*f) * (*g);
                    result.append(newTerm);
                } else {
                    MSG_ABORT("Invalid Gaussian type!");
                }
            } else {
                MSG_ABORT("Invalid Gaussian type!");
            }
        }
    }
    return result;
}

template <int D> GaussExp<D> GaussExp<D>::mult(GaussFunc<D> &g) {
    GaussExp<D> result;
    int nTerms = this->size();
    for (int n = 0; n < nTerms; n++) {
        if (auto *f = dynamic_cast<GaussFunc<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm = *f * g;
            result.append(newTerm);
        } else if (auto *f = dynamic_cast<GaussPoly<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm = *f * g;
            result.append(newTerm);
        } else {
            MSG_ABORT("Invalid Gaussian type!");
        }
    }
    return result;
}
template <int D> GaussExp<D> GaussExp<D>::mult(GaussPoly<D> &g) {
    int nTerms = this->size();
    GaussExp<D> result(nTerms);
    for (int n = 0; n < nTerms; n++) {
        if (auto *f = dynamic_cast<GaussFunc<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm(g * *f);
            result.append(newTerm);
        } else if (auto *f = dynamic_cast<GaussPoly<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm(g * *f);
            result.append(newTerm);
        } else {
            MSG_ABORT("Invalid Gaussian type!");
        }
    }
    return result;
}

template <int D> GaussExp<D> GaussExp<D>::mult(double d) {
    GaussExp<D> prod = *this;
    for (int i = 0; i < this->size(); i++) prod.funcs[i]->multConstInPlace(d);
    return prod;
}

template <int D> void GaussExp<D>::multInPlace(double d) {
    for (int i = 0; i < this->size(); i++) this->funcs[i]->multConstInPlace(d);
}

template <int D> double GaussExp<D>::calcSquareNorm() const {
    /* computing the squares */
    double norm = 0.0;
    for (int i = 0; i < this->size(); i++) {
        double nc = this->funcs[i]->calcSquareNorm();
        norm += nc;
    }
    /* computing the double products */
    for (int i = 0; i < this->size(); i++) {
        GaussExp<D> funcs_i = getFunc(i).asGaussExp(); // Make sure all entries are GaussFunc
        for (int fi = 0; fi < funcs_i.size(); fi++) {
            GaussFunc<D> &func_i = static_cast<GaussFunc<D> &>(funcs_i.getFunc(fi));
            for (int j = i + 1; j < this->size(); j++) {
                GaussExp<D> funcs_j = getFunc(j).asGaussExp(); // Make sure all entries are GaussFunc
                for (int fj = 0; fj < funcs_j.size(); fj++) {
                    GaussFunc<D> &func_j = static_cast<GaussFunc<D> &>(funcs_j.getFunc(fj));
                    double overlap = func_i.calcOverlap(func_j);
                    norm += 2.0 * overlap;
                }
            }
        }
    }
    return norm;
}

template <int D> void GaussExp<D>::normalize() {
    double norm = std::sqrt(this->calcSquareNorm());
    for (int i = 0; i < this->size(); i++) {
        double coef = this->funcs[i]->getCoef();
        this->funcs[i]->setCoef(coef / norm);
    }
}

template <int D> void GaussExp<D>::calcScreening(double nStdDev) {
    screening = nStdDev;
    for (int i = 0; i < this->size(); i++) { this->funcs[i]->calcScreening(nStdDev); }
}

template <int D> void GaussExp<D>::setScreen(bool screen) {
    if (screen) {
        this->screening = std::abs(this->screening);
    } else {
        this->screening = -std::abs(this->screening);
    }
    for (int i = 0; i < this->size(); i++) { this->funcs[i]->setScreen(screen); }
}

// Calculate the scaling and wavelet coefs of all the children, and do the
// outer product to make the nD-scaling coefs. Since a Gaussian expansion
// is not separable, we have to do the projection term by term.
/*
template<int D>
void GaussExp<D>::calcWaveletCoefs(MWNode<D> &node) {
    static const int tDim = 1 << D;
    const ScalingBasis &sf = node.getMWTree().getScalingFunctions();
    MatrixXd &scaling = node.getMWTree().getTmpScalingCoefs();
    VectorXd &tmpvec = node.getMWTree().getTmpScalingVector();
    int kp1 = node.getKp1();
    int kp1_d = node.getKp1_d();
    int inpos = kp1_d - kp1;
    int scale = node.getNodeIndex().scale() + 1;
    node.allocCoefs();
    for (int child = 0; child < tDim; child++) {
        int l[D];
        node.calcChildTranslation(child, l);
        for (int n = 0; n < this->size(); n++) {
            if (this->getFunc(n).checkScreen(scale, l)) {
                continue;
            }
            sf.calcScalingCoefs(this->getFunc(n), scale, l, scaling);
            tmpvec.segment(inpos, kp1) = scaling.col(0);
            math_utils::tensorExpandCoefs(D, 0, kp1, kp1_d, scaling, tmpvec);
            node.getCoefs().segment(child * kp1_d, kp1_d) += tmpvec;
        }
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}
*/

template <int D> void GaussExp<D>::setDefaultScreening(double screen) {
    if (screen < 0) { MSG_ERROR("Screening constant cannot be negative!"); }
    defaultScreening = screen;
}

template <int D> std::ostream &GaussExp<D>::print(std::ostream &o) const {
    o << "Gaussian expansion: " << size() << " terms" << std::endl;
    for (int i = 0; i < size(); i++) {
        o << "Term" << std::setw(3) << i << " :" << std::endl;
        o << getFunc(i) << std::endl << std::endl;
    }
    return o;
}

/** @returns Coulomb repulsion energy between all pairs in GaussExp, including self-interaction
 *
 *  @note Each Gaussian must be normalized to unit charge
 *  \f$ c = (\alpha/\pi)^{D/2} \f$ for this to be correct!
 */
template <int D> double GaussExp<D>::calcCoulombEnergy() const {
    NOT_IMPLEMENTED_ABORT
}

template <> double GaussExp<3>::calcCoulombEnergy() const {
    double energy = 0.0;
    for (int i = 0; i < this->size(); i++) {
        GaussExp<3> funcs_i = getFunc(i).asGaussExp(); // Make sure all entries are GaussFunc
        for (int fi = 0; fi < funcs_i.size(); fi++) {
            GaussFunc<3> &func_i = static_cast<GaussFunc<3> &>(funcs_i.getFunc(fi));
            for (int j = i; j < this->size(); j++) {
                GaussExp<3> funcs_j = getFunc(j).asGaussExp(); // Make sure all entries are GaussFunc
                for (int fj = 0; fj < funcs_j.size(); fj++) {
                    GaussFunc<3> &func_j = static_cast<GaussFunc<3> &>(funcs_j.getFunc(fj));
                    double c = 2.0;
                    if (i == j) c = 1.0;
                    energy += c * func_i.calcCoulombEnergy(func_j);
                }
            }
        }
    }
    return energy;
}

template <int D> GaussExp<D> GaussExp<D>::periodify(const std::array<double, D> &period, double nStdDev) const {
    GaussExp<D> out_exp;
    for (const auto &gauss : *this) {
        auto periodic_gauss = gauss->periodify(period, nStdDev);
        out_exp.append(periodic_gauss);
    }
    return out_exp;
}

template class GaussExp<1>;
template class GaussExp<2>;
template class GaussExp<3>;

} // namespace mrcpp
