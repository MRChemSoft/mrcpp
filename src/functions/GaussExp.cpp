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

template <int D, typename T> double GaussExp<D, T>::defaultScreening = 10.0;

template <int D, typename T> GaussExp<D, T>::GaussExp(int nTerms, double prec) {
    for (int i = 0; i < nTerms; i++) { this->funcs.push_back(nullptr); }
}

template <int D, typename T> GaussExp<D, T>::GaussExp(const GaussExp<D, T> &gexp) {
    screening = gexp.screening;
    for (unsigned int i = 0; i < gexp.size(); i++) {
        Gaussian<D, T> *gauss = gexp.funcs[i]->copy();
        this->funcs.push_back(gauss);
    }
}

template <int D, typename T> GaussExp<D, T>::~GaussExp() {
    for (int i = 0; i < size(); i++) {
        if (this->funcs[i] != nullptr) {
            delete this->funcs[i];
            this->funcs[i] = nullptr;
        }
    }
}

template <int D, typename T> GaussExp<D, T> &GaussExp<D, T>::operator=(const GaussExp<D, T> &gexp) {
    if (&gexp == this) return *this;
    // screening = gexp.screening;
    this->funcs.clear();
    for (unsigned int i = 0; i < gexp.size(); i++) {
        if (gexp.funcs[i] == nullptr) {
            this->funcs.push_back(nullptr);
        } else {
            Gaussian<D, T> *gauss = gexp.getFunc(i).copy();
            this->funcs.push_back(gauss);
        }
    }
    return *this;
}

template <int D, typename T> T GaussExp<D, T>::evalf(const Coord<D> &r) const {
    T val = 0.0;
    for (int i = 0; i < this->size(); i++) { val += this->getFunc(i).evalf(r); }
    return val;
}

template <int D, typename T> bool GaussExp<D, T>::isVisibleAtScale(int scale, int nPts) const {
    for (unsigned int i = 0; i < this->size(); i++) {
        if (not this->getFunc(i).isVisibleAtScale(scale, nPts)) { return false; }
    }
    return true;
}

template <int D, typename T> bool GaussExp<D, T>::isZeroOnInterval(const double *lb, const double *ub) const {
    for (unsigned int i = 0; i < this->size(); i++) {
        if (not this->getFunc(i).isZeroOnInterval(lb, ub)) { return false; }
    }
    return true;
}

template <int D, typename T> void GaussExp<D, T>::setFunc(int i, const GaussPoly<D, T> &g, double c) {
    if (i < 0 or i > (this->size() - 1)) {
        MSG_ERROR("Index out of bounds!");
        return;
    }
    if (this->funcs[i] != nullptr) { delete this->funcs[i]; }
    this->funcs[i] = new GaussPoly<D, T>(g);
    double coef = this->funcs[i]->getCoef();
    this->funcs[i]->setCoef(c * coef);
}

template <int D, typename T> void GaussExp<D, T>::setFunc(int i, const GaussFunc<D, T> &g, double c) {
    if (i < 0 or i > (this->size() - 1)) {
        MSG_ERROR("Index out of bounds!");
        return;
    }
    if (this->funcs[i] != nullptr) { delete this->funcs[i]; }
    this->funcs[i] = new GaussFunc<D, T>(g);
    double coef = this->funcs[i]->getCoef();
    this->funcs[i]->setCoef(c * coef);
}

template <int D, typename T> void GaussExp<D, T>::append(const Gaussian<D, T> &g) {
    Gaussian<D, T> *gp = g.copy();
    this->funcs.push_back(gp);
}

template <int D, typename T> void GaussExp<D, T>::append(const GaussExp<D, T> &g) {
    for (int i = 0; i < g.size(); i++) {
        Gaussian<D, T> *gp = g.getFunc(i).copy();
        this->funcs.push_back(gp);
    }
}

template <int D, typename T> GaussExp<D, T> GaussExp<D, T>::differentiate(int dir) const {
    assert(dir >= 0 and dir < D);
    GaussExp<D, T> result;
    for (int i = 0; i < this->size(); i++) result.append(this->getFunc(i).differentiate(dir));
    return result;
}

template <int D, typename T> GaussExp<D, T> GaussExp<D, T>::add(GaussExp<D, T> &g) {
    int nsum = this->size() + g.size();
    GaussExp<D, T> sum = GaussExp<D, T>(nsum);

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

template <int D, typename T> GaussExp<D, T> GaussExp<D, T>::add(Gaussian<D, T> &g) {
    int nsum = this->size() + 1;
    GaussExp<D, T> sum = GaussExp<D, T>(nsum);
    for (int n = 0; n < this->size(); n++) { sum.funcs[n] = this->getFunc(n).copy(); }
    sum.funcs[this->size()] = g.copy();
    return sum;
}

template <int D, typename T> GaussExp<D, T> GaussExp<D, T>::mult(GaussExp<D, T> &gexp) {
    GaussExp<D, T> result;
    for (int i = 0; i < this->size(); i++) {
        for (int j = 0; j < gexp.size(); j++) {
            if (auto *f = dynamic_cast<GaussFunc<D, T> *>(this->funcs[i])) {
                if (auto *g = dynamic_cast<GaussFunc<D, T> *>(gexp.funcs[j])) {
                    GaussPoly<D, T> newTerm = (*g) * (*f);
                    result.append(newTerm);
                } else if (auto *g = dynamic_cast<GaussPoly<D, T> *>(gexp.funcs[j])) {
                    GaussPoly<D, T> newTerm = (*g) * (*f);
                    result.append(newTerm);
                } else {
                    MSG_ABORT("Invalid Gaussian type!");
                }
            } else if (auto *f = dynamic_cast<GaussPoly<D, T> *>(this->funcs[i])) {
                if (auto *g = dynamic_cast<GaussFunc<D, T> *>(gexp.funcs[j])) {
                    GaussPoly<D, T> newTerm = (*f) * (*g);
                    result.append(newTerm);
                } else if (auto *g = dynamic_cast<GaussPoly<D, T> *>(gexp.funcs[j])) {
                    GaussPoly<D, T> newTerm = (*f) * (*g);
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

template <int D, typename T> GaussExp<D, T> GaussExp<D, T>::mult(GaussFunc<D, T> &g) {
    GaussExp<D, T> result;
    int nTerms = this->size();
    for (int n = 0; n < nTerms; n++) {
        if (auto *f = dynamic_cast<GaussFunc<D, T> *>(this->funcs[n])) {
            GaussPoly<D, T> newTerm = *f * g;
            result.append(newTerm);
        } else if (auto *f = dynamic_cast<GaussPoly<D, T> *>(this->funcs[n])) {
            GaussPoly<D, T> newTerm = *f * g;
            result.append(newTerm);
        } else {
            MSG_ABORT("Invalid Gaussian type!");
        }
    }
    return result;
}
template <int D, typename T> GaussExp<D, T> GaussExp<D, T>::mult(GaussPoly<D, T> &g) {
    int nTerms = this->size();
    GaussExp<D, T> result(nTerms);
    for (int n = 0; n < nTerms; n++) {
        if (auto *f = dynamic_cast<GaussFunc<D, T> *>(this->funcs[n])) {
            GaussPoly<D, T> newTerm(g * *f);
            result.append(newTerm);
        } else if (auto *f = dynamic_cast<GaussPoly<D, T> *>(this->funcs[n])) {
            GaussPoly<D, T> newTerm(g * *f);
            result.append(newTerm);
        } else {
            MSG_ABORT("Invalid Gaussian type!");
        }
    }
    return result;
}

template <int D, typename T> GaussExp<D, T> GaussExp<D, T>::mult(double d) {
    GaussExp<D, T> prod = *this;
    for (int i = 0; i < this->size(); i++) prod.funcs[i]->multConstInPlace(d);
    return prod;
}

template <int D, typename T> void GaussExp<D, T>::multInPlace(double d) {
    for (int i = 0; i < this->size(); i++) this->funcs[i]->multConstInPlace(d);
}

template <int D, typename T> double GaussExp<D, T>::calcSquareNorm() const {
    /* computing the squares */
    double norm = 0.0;
    for (int i = 0; i < this->size(); i++) {
        double nc = this->funcs[i]->calcSquareNorm();
        norm += nc;
    }
    /* computing the double products */
    for (int i = 0; i < this->size(); i++) {
        GaussExp<D, T> funcs_i = getFunc(i).asGaussExp(); // Make sure all entries are GaussFunc
        for (int fi = 0; fi < funcs_i.size(); fi++) {
            GaussFunc<D, T> &func_i = static_cast<GaussFunc<D, T> &>(funcs_i.getFunc(fi));
            for (int j = i + 1; j < this->size(); j++) {
                GaussExp<D, T> funcs_j = getFunc(j).asGaussExp(); // Make sure all entries are GaussFunc
                for (int fj = 0; fj < funcs_j.size(); fj++) {
                    GaussFunc<D, T> &func_j = static_cast<GaussFunc<D, T> &>(funcs_j.getFunc(fj));
                    double overlap = func_i.calcOverlap(func_j);
                    norm += 2.0 * overlap;
                }
            }
        }
    }
    return norm;
}

template <int D, typename T> void GaussExp<D, T>::normalize() {
    double norm = std::sqrt(this->calcSquareNorm());
    for (int i = 0; i < this->size(); i++) {
        double coef = this->funcs[i]->getCoef();
        this->funcs[i]->setCoef(coef / norm);
    }
}

template <int D, typename T> void GaussExp<D, T>::calcScreening(double nStdDev) {
    screening = nStdDev;
    for (int i = 0; i < this->size(); i++) { this->funcs[i]->calcScreening(nStdDev); }
}

template <int D, typename T> void GaussExp<D, T>::setScreen(bool screen) {
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
void GaussExp<D, T>::calcWaveletCoefs(MWNode<D, T> &node) {
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

template <int D, typename T> void GaussExp<D, T>::setDefaultScreening(double screen) {
    if (screen < 0) { MSG_ERROR("Screening constant cannot be negative!"); }
    defaultScreening = screen;
}

template <int D, typename T> std::ostream &GaussExp<D, T>::print(std::ostream &o) const {
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
template <int D, typename T> double GaussExp<D, T>::calcCoulombEnergy() const {
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

template <int D, typename T> GaussExp<D, T> GaussExp<D, T>::periodify(const std::array<double, D> &period, double nStdDev) const {
    GaussExp<D, T> out_exp;
    for (const auto &gauss : *this) {
        auto periodic_gauss = gauss->periodify(period, nStdDev);
        out_exp.append(periodic_gauss);
    }
    return out_exp;
}

template class GaussExp<1, double>;
template class GaussExp<2, double>;
template class GaussExp<3, double>;

template class GaussExp<1, ComplexDouble>;
template class GaussExp<2, ComplexDouble>;
template class GaussExp<3, ComplexDouble>;

} // namespace mrcpp
