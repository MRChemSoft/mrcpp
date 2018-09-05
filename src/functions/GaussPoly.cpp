/**
 *
 *
 * \date May 26, 2010
 * \author Stig Rune Jensen
 *		   CTCC, University of Tromsø
 *
 */

#include "vector"

#include "GaussPoly.h"
#include "GaussFunc.h"
#include "GaussExp.h"
#include "utils/Printer.h"

using namespace std;
using namespace Eigen;

namespace mrcpp {

template<int D>
GaussPoly<D>::GaussPoly(double alpha, double coef, const double pos[D],
    const int pow[D]) : Gaussian<D>(alpha, coef, pos, pow) {
    for (int d = 0; d < D; d++) {
        if (pow != 0) {
            this->poly[d] = new Polynomial(this->power[d]);
            //this->poly[d]->unsetBounds();
        } else {
            this->poly[d] = 0;
        }
    }
}

template<int D>
GaussPoly<D>::GaussPoly(const GaussPoly<D> &gp) : Gaussian<D>(gp) {
    for (int d = 0; d < D; d++) {
        poly[d] = new Polynomial(gp.getPoly(d));
    }
}

template<int D>
GaussPoly<D>::GaussPoly(const GaussFunc<D> &gf) : Gaussian<D>(gf) {
    for (int d = 0; d < D; d++) {
        int order = this->getPower(d);
        poly[d] = new Polynomial(order);
        VectorXd coefs = VectorXd::Zero(order + 1);
        coefs[order] = 1.0;
        poly[d]->setCoefs(coefs);
        //poly[d]->unsetBounds();
    }
}

template<int D>
GaussPoly<D>::~GaussPoly() {
    for (int i = 0; i < D; i++) {
        delete poly[i];
    }
}

template<int D>
Gaussian<D> *GaussPoly<D>::copy() const {
    GaussPoly<D> *gauss = new GaussPoly<D>(*this);
    return gauss;
}

template<int D>
double GaussPoly<D>::calcOverlap(GaussPoly<D> &b) {
    GaussExp<D> gExp(*this);
    GaussExp<D> fExp(b);

    double overlap = 0.0;
    for (int i = 0; i < fExp.size(); i++) {
        GaussFunc<D> &fFunc = static_cast<GaussFunc<D> &>(fExp.getFunc(i));
        for (int j = 0; j < gExp.size(); j++) {
            overlap += gExp.getFunc(j).calcOverlap(fFunc);
        }
    }
    return overlap;
}

template<int D>
double GaussPoly<D>::calcOverlap(GaussFunc<D> &b) {
    return b.calcOverlap(*this);
}

template<int D>
double GaussPoly<D>::calcSquareNorm() {
    this->squareNorm = this->calcOverlap(*this);
    return this->squareNorm;
}

template<int D>
double GaussPoly<D>::evalf(const double *r) const {
    if (this->getScreen()) {
        for (int d = 0; d < D; d++) {
            if ((r[d] < this->A[d]) or (r[d] > this->B[d])) {
                return 0.0;
            }
        }
    }
    double q2 = 0.0, p2 = 1.0;
    for (int d = 0; d < D; d++) {
        //assert(this->poly[d]->getCheckBounds() == false);
        double q = r[d] - this->pos[d];
        q2 += q * q;
        p2 *= poly[d]->evalf(r[d] - this->pos[d]);
    }
    return this->coef * p2 * exp(-this->alpha * q2);
}

/** NOTE!
 *	This function evaluation will give the first dimension the full coef
 *	amplitude, leaving all other directions with amplitude 1.0. This is to
 *	avoid expensive d-root evaluation when distributing the amplitude
 *	equally to all dimensions.
 */
template<int D>
double GaussPoly<D>::evalf(const double r, int d) const {
    if (this->getScreen()) {
        if ((r < this->A[d]) or (r > this->B[d])) {
            return 0.0;
        }
    }
    //assert(this->poly[d]->getCheckBounds() == false);
    double q2 = 0.0, p2 = 1.0;
    double q = (r - this->pos[d]);
    q2 += q * q;
    p2 *= poly[d]->evalf(q);
    if (d == 0) {
        p2 *= this->coef;
    }
    return p2 * exp(-this->alpha * q2);
}

template<int D>
GaussPoly<D> GaussPoly<D>::differentiate(int dir) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GaussPoly<D>::multInPlace(const GaussPoly<D> &rhs) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GaussPoly<D>::fillCoefPowVector(vector<double> &coefs, vector<int *> &power,
        int pow[D], int dir) const {
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

template<int D>
GaussPoly<D> GaussPoly<D>::mult(const GaussPoly<D> &rhs) {
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

template<int D>
GaussPoly<D> GaussPoly<D>::mult(double c) {
    GaussPoly<D> g = *this;
    g.coef *= c;
    return g;
}

template<int D>
void GaussPoly<D>::setPower(int d, int pow) {
    if (poly[d] != 0) {
        delete poly[d];
    }
    poly[d] = new Polynomial(pow);
    this->squareNorm = -1.0;
}

template<int D>
void GaussPoly<D>::setPower(const int pow[D]) {
    for (int d = 0; d < D; d++) {
        if (poly[d] != 0) {
            delete poly[d];
        }
        poly[d] = new Polynomial(pow[d]);
    }
    this->squareNorm = -1.0;
}

template<int D>
void GaussPoly<D>::setPoly(int d, Polynomial &poly) {
    if (this->poly[d] != 0) {
        delete this->poly[d];
    }
    this->poly[d] = new Polynomial(poly);
    //this->poly[d]->unsetBounds();
    this->power[d] = poly.getOrder();
}

template<int D>
std::ostream& GaussPoly<D>::print(std::ostream &o) const {
    o << "Exp: " << this->getExp() << std::endl;
    o << "Coef: " << this->getCoef() << std::endl;
    o << "Pos:   ";
    for (int i = 0; i < D; i++) {
        o << this->getPos()[i] << " ";
    }
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
