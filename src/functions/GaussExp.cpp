/**
 */

// TODO checks on exponents to discard elements
// TODO sanity checks

#include <iostream>
#include <cstdlib>

#include "GaussExp.h"
#include "Gaussian.h"
#include "GaussFunc.h"
#include "GaussPoly.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

using namespace std;
using namespace Eigen;

namespace mrcpp {

template<int D> double GaussExp<D>::defaultScreening = 10.0;

template<int D>
GaussExp<D>::GaussExp(int nTerms, double prec) : screening(0.0), squareNorm(-1.0) {
    for (int i = 0; i < nTerms; i++) {
        this->funcs.push_back(0);
    }
}

template<int D>
GaussExp<D>::GaussExp(const GaussExp<D> &gexp) {
    this->squareNorm = gexp.squareNorm;
    screening = gexp.screening;
    for (unsigned int i = 0; i < gexp.size(); i++) {
        Gaussian<D> *gauss = gexp.funcs[i]->copy();
        this->funcs.push_back(gauss);
    }
}

template<int D>
GaussExp<D>::GaussExp(const GaussPoly<D> &gPoly) : screening(0.0), squareNorm(-1.0) {
    int pow[D];
    double pos[D];
    double coef;
    double alpha = gPoly.getExp();

    int nTerms = 1;
    for (int d = 0; d < D; d++) {
        nTerms *= (gPoly.getPower(d) + 1);
        pos[d] = gPoly.getPos()[d];
    }

    vector<double> coefs;
    vector<int *> power;

    gPoly.fillCoefPowVector(coefs, power, pow, D);

    for (int i = 0; i < nTerms; i++) {
        coef = coefs[i];
        for (int d = 0; d < D; d++) {
            pow[d] = power[i][d];
        }
        if (coef != 0.0) {
            GaussFunc<D> gFunc(alpha, coef, pos, pow);
            this->append(gFunc);
        }
    }
    for (int i = 0; i < power.size(); i++) {
        delete[] power[i];
    }
}

template<int D>
GaussExp<D>::~GaussExp() {
    for (int i = 0; i < size(); i++) {
        if (this->funcs[i] != 0) {
            delete this->funcs[i];
            this->funcs[i] = 0;
        }
    }
}

template<int D>
GaussExp<D> &GaussExp<D>::operator=(const GaussExp<D> &gexp) {
    if (&gexp == this)
        return *this;
    squareNorm = gexp.squareNorm;
    //screening = gexp.screening;
    this->funcs.clear();
    for (unsigned int i = 0; i < gexp.size(); i++) {
        if (gexp.funcs[i] == 0) {
            this->funcs.push_back(0);
        } else {
            Gaussian<D> *gauss = gexp.getFunc(i).copy();
            this->funcs.push_back(gauss);
        }
    }
    return *this;
}

template<int D>
double GaussExp<D>::evalf(const double *r) const {
    double val = 0.0;
    for (int i = 0; i < this->size(); i++) {
        val += this->getFunc(i).evalf(r);
    }
    return val;
}

template<int D>
bool GaussExp<D>::isVisibleAtScale(int scale, int nPts) const {
    for (unsigned int i = 0; i < this->size(); i++) {
        if (not this->getFunc(i).isVisibleAtScale(scale, nPts)) {
            return false;
        }
    }
    return true;
}

template<int D>
bool GaussExp<D>::isZeroOnInterval(const double *lb, const double *ub) const {
    for (unsigned int i = 0; i < this->size(); i++) {
        if (not this->getFunc(i).isZeroOnInterval(lb, ub)) {
            return false;
        }
    }
    return true;
}

template<int D>
void GaussExp<D>::setFunc(int i, const GaussPoly<D> &g, double c) {
    if (i < 0 or i > (this->size() - 1)) {
        MSG_ERROR("Index out of bounds!");
        return;
    }
    if (this->funcs[i] != 0) {
        delete this->funcs[i];
    }
    this->funcs[i] = new GaussPoly<D> (g);
    double coef = this->funcs[i]->getCoef();
    this->funcs[i]->setCoef(c*coef);
}

template<int D>
void GaussExp<D>::setFunc(int i, const GaussFunc<D> &g, double c) {
    if (i < 0 or i > (this->size() - 1)) {
        MSG_ERROR("Index out of bounds!");
        return;
    }
    if (this->funcs[i] != 0) {
        delete this->funcs[i];
    }
    this->funcs[i] = new GaussFunc<D> (g);
    double coef = this->funcs[i]->getCoef();
    this->funcs[i]->setCoef(c*coef);
}

template<int D>
void GaussExp<D>::append(const Gaussian<D> &g) {
    Gaussian<D> *gp = g.copy();
    this->funcs.push_back(gp);
    this->squareNorm = -1.0;
}

template<int D>
void GaussExp<D>::append(const GaussExp<D> &g) {
    for (int i = 0; i < g.size(); i++) {
        Gaussian<D> *gp = g.getFunc(i).copy();
        this->funcs.push_back(gp);
    }
    this->squareNorm = -1.0;
}

template<int D>
GaussExp<D> GaussExp<D>::differentiate(int dir) {
    assert(dir >= 0 and dir < D);
    GaussExp<D> result;
    for (int i = 0; i < this->size(); i++) {
        result.append(this->getFunc(i).differentiate(dir));
    }
    return result;
}

template<int D>
GaussExp<D> GaussExp<D>::add(GaussExp<D> & g) {
    int nsum = this->size() + g.size();
    GaussExp<D> sum = GaussExp<D> (nsum);

    int n = 0;
    for (int i = 0; i < this->size(); i++) {
        sum.funcs[n] = this->funcs[i]->copy();
        n++;
    }
    for (int i = 0; i < g.size(); i++) {
        sum.funcs[n] = g.funcs[i]->copy();
        n++;
    }

    sum.calcSquareNorm();

    return sum;
}

template<int D>
GaussExp<D> GaussExp<D>::add(Gaussian<D> &g) {
    int nsum = this->size() + 1;
    GaussExp<D> sum = GaussExp<D> (nsum);
    for (int n = 0; n < this->size(); n++) {
        sum.funcs[n] = this->getFunc(n).copy();
    }
    sum.funcs[this->size()] = g.copy();
    return sum;
}

template<int D>
GaussExp<D> GaussExp<D>::mult(GaussExp<D> &gexp) {
    GaussExp<D> result;
    for (int i = 0; i < this->size(); i++) {
        for (int j = 0; j < gexp.size(); j++) {
            if (GaussFunc<D> *f =
                                dynamic_cast<GaussFunc<D> *>(this->funcs[i])) {
                if (GaussFunc<D> *g =
                                dynamic_cast<GaussFunc<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*g) * (*f);
                    result.append(newTerm);
                } else if (GaussPoly<D> *g =
                                dynamic_cast<GaussPoly<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*g) * (*f);
                    result.append(newTerm);
                } else {
                    MSG_FATAL("Invalid Gaussian type!");
                }
            } else if (GaussPoly<D> *f =
                                dynamic_cast<GaussPoly<D> *>(this->funcs[i])) {
                if (GaussFunc<D> *g =
                                dynamic_cast<GaussFunc<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*f) * (*g);
                    result.append(newTerm);
                } else if (GaussPoly<D> *g =
                                dynamic_cast<GaussPoly<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*f) * (*g);
                    result.append(newTerm);
                } else {
                    MSG_FATAL("Invalid Gaussian type!");
                }
            } else {
                MSG_FATAL("Invalid Gaussian type!");
            }
        }
    }
    return result;
}

template<int D>
GaussExp<D> GaussExp<D>::mult(GaussFunc<D> &g) {
    GaussExp<D> result;
    int nTerms = this->size();
    for (int n = 0; n < nTerms; n++) {
        if (GaussFunc<D> *f = dynamic_cast<GaussFunc<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm = *f * g;
            result.append(newTerm);
        } else if (GaussPoly<D> *f = dynamic_cast<GaussPoly<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm = *f * g;
            result.append(newTerm);
        } else {
            MSG_FATAL("Invalid Gaussian type!");
        }
    }
    return result;
}
template<int D>
GaussExp<D> GaussExp<D>::mult(GaussPoly<D> &g) {
    int nTerms = this->size();
    GaussExp<D> result(nTerms);
    for (int n = 0; n < nTerms; n++) {
        if (GaussFunc<D> *f = dynamic_cast<GaussFunc<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm(g * *f);
            result.append(newTerm);
        } else if (GaussPoly<D> *f = dynamic_cast<GaussPoly<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm(g * *f);
            result.append(newTerm);
        } else {
            MSG_FATAL("Invalid Gaussian type!");
        }
    }
    return result;
}

template<int D>
GaussExp<D> GaussExp<D>::mult(double d) {
    GaussExp<D> prod = *this;

    for (int i = 0; i < this->size(); i++) {
        prod.funcs[i]->multConstInPlace(d);
    }
    prod.calcSquareNorm();
    return prod;
}

template<int D>
void GaussExp<D>::multInPlace(double d) {
    for (int i = 0; i < this->size(); i++) {
        this->funcs[i]->multConstInPlace(d);
    }
//	this->squareNorm = -1.0;
    this->calcSquareNorm();
}

/**  compute the norm of a multidimensional gaussian expansion
 */
template<int D>
double GaussExp<D>::calcSquareNorm() {

    double overlap, norm = 0.0;

    /* computing the squares */
    for (int i = 0; i < this->size(); i++) {
        double nc = this->funcs[i]->getSquareNorm();
        norm += nc;
    }
    /* computing the double products */
    for (int i = 0; i < this->size(); i++) {
        for (int j = i + 1; j < this->size(); j++) {
            if (GaussFunc<D> *f = dynamic_cast<GaussFunc<D> *>(this->funcs[j])) {
                overlap = this->funcs[i]->calcOverlap(*f);
            } else if (GaussPoly<D> *f =
                    dynamic_cast<GaussPoly<D> *>(this->funcs[j])) {
                overlap = this->funcs[i]->calcOverlap(*f);
            } else {
                MSG_FATAL("Invald argument");
            }
            norm += 2.0 * overlap;
        }
    }
    this->squareNorm = norm;
    return norm;
}

template<int D>
void GaussExp<D>::normalize() {
    double norm = sqrt(this->getSquareNorm());
    double coef;
    for (int i = 0; i < this->size(); i++) {
        coef = this->funcs[i]->getCoef();
        this->funcs[i]->setCoef(coef / norm);
    }
    calcSquareNorm();
}

/** calculate screening radius on the individual 1d terms of the gauss
 * expansion
 *
 *	When projecting a large gaussian expansion on a given node we want to
 *	discard the terms in the expansion that are sufficiently far from the
 *	node. These terms will not contribute to the function on that node, but
 *	nevertheless it will take a considerable amount of time to evaluate them.
 *
 *	By applying this routine a screening of far away gaussians is enabled, and
 *	the projection becomes more efficient.
 */
template<int D>
void GaussExp<D>::calcScreening(double nStdDev) {
    screening = nStdDev;
    for (int i = 0; i < this->size(); i++) {
        this->funcs[i]->calcScreening(nStdDev);
    }
}

template<int D>
void GaussExp<D>::setScreen(bool screen) {
    if (screen) {
        this->screening = fabs(this->screening);
    } else {
        this->screening = -fabs(this->screening);
    }
    for (int i = 0; i < this->size(); i++) {
        this->funcs[i]->setScreen(screen);
    }
}


/** Calculate the scaling and wavelet coefs of all the children, and do the
 * outer product to make the nD-scaling coefs. Since a Gaussian expansion
 * is not separable, we have to do the projection term by term. */
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

// Specialized for D=3 below
template<int D>
double GaussExp<D>::calcCoulombEnergy() {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void GaussExp<D>::setDefaultScreening(double screen) {
    if (screen < 0) {
        MSG_ERROR("Screening constant cannot be negative!");
    }
    defaultScreening = screen;
}

template<int D>
std::ostream& GaussExp<D>::print(std::ostream &o) const {
    o << "Gaussian Expansion: " << size() << " terms" << std::endl;
    for (int i = 0; i < size(); i++) {
        o << "Term " << i << ":" << std::endl;
        o << getFunc(i) << std::endl << std::endl;
    }
    return o;
}

template<>
double GaussExp<3>::calcCoulombEnergy() {
    double energy = 0.0;
    for (int i = 0; i < size(); i++) {
        if (GaussFunc<3> *gauss_i = dynamic_cast<GaussFunc<3> *>(&getFunc(i))) {
            for (int j = i; j < size(); j++) {
                if (GaussFunc<3> *gauss_j = dynamic_cast<GaussFunc<3> *>(&getFunc(j))) {
                    double c = 2.0;
                    if (i == j) c = 1.0;
                    energy += c*gauss_i->calcCoulombEnergy(*gauss_j);
                } else {
                    MSG_ERROR("Can only calculate energy for GaussFunc");
                }
            }
        } else {
            MSG_ERROR("Can only calculate energy for GaussFunc");
        }
    }
    return energy;
}

template class GaussExp<1>;
template class GaussExp<2>;
template class GaussExp<3>;

} // namespace mrcpp
