/*
 *
 *
 *  \date Jul 5, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */
#include <cmath>

#include "GaussExp.h"
#include "GaussFunc.h"
#include "GaussPoly.h"
#include "Polynomial.h"
#include "BoysFunction.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

template<int D>
Gaussian<D> *GaussFunc<D>::copy() const{
    GaussFunc<D> *gauss = new GaussFunc<D>(*this);
    return gauss;
}

template<int D>
double GaussFunc<D>::evalf(const double *r) const {
    if (this->getScreen()) {
        for (int d = 0; d < D; d++) {
            if (r[d] < this->A[d] or r[d] > this->B[d]) {
                return 0.0;
            }
        }
    }
    double q2 = 0.0, p2 = 1.0;
    for (int d = 0; d < D; d++) {
        double q = r[d] - this->pos[d];
        q2 += q * q;
        if (this->power[d] == 0) {
            continue;
        } else if (this->power[d] == 1) {
            p2 *= q;
        } else {
            p2 *= pow(q, this->power[d]);
        }
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
double GaussFunc<D>::evalf(double r, int d) const {
    if (this->getScreen()) {
        if ((r < this->A[d]) or (r > this->B[d])) {
            return 0.0;
        }
    }
    double q;
    double q2, p2;
    q = (r - this->pos[d]);
    q2 = q * q;
    if (this->power[d] == 0) {
        p2 = 1.0;
    } else if (this->power[d] == 1) {
        p2 = q;
    } else {
        p2 = pow(q, this->power[d]);
    }
    double result = p2 * exp(-this->alpha * q2);
    if (d == 0)
        result *= this->coef;
    return result;
}

template<int D>
double GaussFunc<D>::calcSquareNorm() {
    int p, i;
    double sq_norm, a;
    double norm = 1.0;

    for (int n = 0; n < D; n++) {
        a = 2.0 * this->alpha;
        sq_norm = 1.0;
        p = this->power[n];
        if (p > 0) {
            i = 2 * p - 1;
            while (i > 0) {
                sq_norm = i * sq_norm / (2.0 * a);
                i = i - 2;
            }
        }
        a = pi / a;
        sq_norm *= sqrt(a);
        norm *= sq_norm;
    }
    this->squareNorm = norm * this->coef * this->coef;
    return this->squareNorm;
}

template<int D>
GaussPoly<D> GaussFunc<D>::differentiate(int dir) {
    GaussPoly<D> result(*this);
    int oldPow = this->getPower(dir);

    Polynomial newPoly(oldPow + 1);
    newPoly.getCoefs()[oldPow + 1] = -2.0 * this->getExp();
    if (oldPow > 0) {
        newPoly.getCoefs()[oldPow - 1] = oldPow;
    }
    result.setPoly(dir, newPoly);;
    return result;
}

template<int D>
void GaussFunc<D>::multInPlace(const GaussFunc<D> &rhs) {
    GaussFunc<D> &lhs = *this;
    for (int d = 0; d < D; d++) {
        if (lhs.getPos()[d] != rhs.getPos()[d]) {
            MSG_FATAL("Cannot multiply GaussFuncs of different center in-place");
        }
    }
    double newCoef = lhs.getCoef() * rhs.getCoef();
    double newExp = lhs.getExp() + rhs.getExp();
    int newPow[D];
    for (int d = 0; d < D; d++) {
        newPow[d] = lhs.getPower(d) + rhs.getPower(d);
    }
    this->setCoef(newCoef);
    this->setExp(newExp);
    this->setPower(newPow);
//	this->squareNorm = -1.0;
    this->calcSquareNorm();
}

template<int D>
GaussPoly<D> GaussFunc<D>::mult(const GaussFunc<D> &rhs) {
    GaussFunc<D> &lhs = *this;
    GaussPoly<D> result;
    result.multPureGauss(lhs, rhs);
    for (int d = 0; d < D; d++) {
        double newPos = result.getPos()[d];
        Polynomial lhsPoly(newPos - lhs.getPos()[d], lhs.getPower(d));
        Polynomial rhsPoly(newPos - rhs.getPos()[d], rhs.getPower(d));
        Polynomial newPoly = lhsPoly * rhsPoly;
        result.setPoly(d, newPoly);
    }
    result.setCoef(result.getCoef() * lhs.getCoef() * rhs.getCoef());
    return result;
}

template<int D>
GaussFunc<D> GaussFunc<D>::mult(double c) {
    GaussFunc<D> g = *this;
    g.coef *= c;
    return g;
}

template<int D>
double GaussFunc<D>::calcOverlap(GaussPoly<D> &b) {
    GaussExp<D> gExp(b);
    double overlap = 0.0;
    for (int i = 0; i < gExp.size(); i++) {
        GaussFunc<D> &gFunc = static_cast<GaussFunc<D> &>(gExp.getFunc(i));
        overlap += this->calcOverlap(gFunc);
    }
    return overlap;
}

/**  Compute the D-dimensional overlap integral between two
 * gaussian distributions */
template<int D>
double GaussFunc<D>::calcOverlap(GaussFunc<D> &b) {
    double S = 1.0;
    for (int d = 0; d < D; d++) {
        S *= ObaraSaika_ab(this->power[d], b.power[d], this->pos[d], b.pos[d],
                this->alpha, b.alpha);
    }
    S *= this->coef * b.coef;
    return S;
}

template<int D>
double GaussFunc<D>::calcOverlap(GaussFunc<D> &a, GaussFunc<D> &b) {
    double S = 1.0;
    for (int d = 0; d < D; d++) {
        S *= ObaraSaika_ab(a.power[d], b.power[d], a.pos[d], b.pos[d], a.alpha,
                b.alpha);
    }
    S *= a.coef * b.coef;
    return S;
}

/**  Compute the monodimensional overlap integral between two
 gaussian distributions by means of the Obara-Saika recursiive
 scheme

 \f[ S_{ij} = \int_{-\infty}^{+\infty} \,\mathrm{d} x
 (x-x_a)^{p_a}
 (x-x_b)^{p_b}
 e^{-c_a (x-x_a)^2}
 e^{-c_b (x-x_b)^2}\f]

 @param power_a \f$ p_a     \f$
 @param power_b \f$ p_b     \f$
 @param pos_a   \f$ x_a     \f$
 @param pos_b   \f$ x_b     \f$
 @param expo_a  \f$ c_a \f$
 @param expo_b  \f$ c_b \f$

 */
template<int D>
double GaussFunc<D>::ObaraSaika_ab(int power_a, int power_b, double pos_a,
        double pos_b, double expo_a, double expo_b) {
    int i, j, i_l, i_r, n_0j_coeff, n_ij_coeff;
    double expo_p, mu, pos_p, x_ab, x_pa, x_pb, s_00;
    /* The highest angular momentum combination is l=20 for a and b
     * simulatnelusly */
    double s_coeff[64];

    //	if (out_of_bounds(power_a, 0, MAX_GAUSS_POWER) ||
    //		out_of_bounds(power_b, 0, MAX_GAUSS_POWER)
    //		) {
    //		PRINT_FUNC_NAME;
    //		INVALID_ARG_EXIT;
    //	}

    /* initialization of a hell of a lot of coefficients.... */
    expo_p = expo_a + expo_b; /* total exponent */
    mu = expo_a * expo_b / (expo_a + expo_b); /* reduced exponent */
    pos_p = (expo_a * pos_a + expo_b * pos_b) / expo_p; /* center of charge */
    x_ab = pos_a - pos_b; /* X_{AB} */
    x_pa = pos_p - pos_a; /* X_{PA} */
    x_pb = pos_p - pos_b; /* X_{PB} */
    s_00 = pi / expo_p;
    s_00 = sqrt(s_00) * exp(-mu * x_ab * x_ab); /* overlap of two spherical gaussians */
    n_0j_coeff = 1 + power_b; /* n. of 0j coefficients needed */
    n_ij_coeff = 2 * power_a; /* n. of ij coefficients needed (i > 0) */

    /* we add 3 coeffs. to avoid a hell of a lot of if statements */
    /*    n_tot_coeff = n_0j_coeff + n_ij_coeff + 3;	*/
    /*    s_coeff = (double *) calloc(n_tot_coeff, sizeof(double));*/

    /* generate first two coefficients */
    s_coeff[0] = s_00;
    s_coeff[1] = x_pb * s_00;
    j = 1;
    /* generate the rest of the first row */
    while (j < power_b) {
        s_coeff[j + 1] = x_pb * s_coeff[j] + j * s_coeff[j - 1]
                / (2.0 * expo_p);
        j++;
    }
    /* generate the first two coefficients with i > 0 */
    s_coeff[j + 1] = s_coeff[j] - x_ab * s_coeff[j - 1];
    s_coeff[j + 2] = x_pa * s_coeff[j] + j * s_coeff[j - 1] / (2.0 * expo_p);
    i = 1;
    /* generate the remaining coefficients with i > 0 */
    while (i < power_a) {
        i_l = j + 2 * i + 1;
        i_r = j + 2 * i + 2;
        s_coeff[i_l] = s_coeff[i_l - 1] - x_ab * s_coeff[i_l - 2];
        s_coeff[i_r] = x_pa * s_coeff[i_r - 2] + (j * s_coeff[i_r - 3] + i
                * s_coeff[i_r - 4]) / (2.0 * expo_p);
        i++;
    }

    /*    free(s_coeff);*/
    return s_coeff[power_b + 2 * power_a];
}

// Specialized for D=3 below
template<int D>
double GaussFunc<D>::calcCoulombEnergy(GaussFunc<D> &gf) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
std::ostream& GaussFunc<D>::print(std::ostream &o) const {
    o << "Exp:   " << this->getExp() << std::endl;
    o << "Coef:  "<< this->getCoef() << std::endl;
    o << "Pos:   ";
    for (int i = 0; i < D; i++) {
        o << this->getPos()[i] << " ";
    }
    o << std::endl;
    o << "Power: ";
    for (int i = 0; i < D; i++) {
        o << this->getPower(i) << " ";
    }
    return o;
}

/** NOTE: Gaussians must be normalized to unit charge coef = (alpha/pi)^(3/2)
 * for this to be correct!
 */
template<>
double GaussFunc<3>::calcCoulombEnergy(GaussFunc<3> &gf) {
    double p = this->getExp();
    double q = gf.getExp();
    double alpha = p*q/(p+q);

    const double *Rp = this->getPos();
    const double *Rq = gf.getPos();

    double Rx = Rp[0] - Rq[0];
    double Ry = Rp[1] - Rq[1];
    double Rz = Rp[2] - Rq[2];

    double Rpq_2 = Rx*Rx + Ry*Ry + Rz*Rz;

    BoysFunction boys(0);

    double boysArg = alpha*Rpq_2;
    double boysFac = boys.evalf(&boysArg);

    return sqrt(4.0*alpha/pi)*boysFac;
}

template class GaussFunc<1>;
template class GaussFunc<2>;
template class GaussFunc<3>;

} // namespace mrcpp
