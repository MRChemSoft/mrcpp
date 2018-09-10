/*
 *
 *  \date Jul 5, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "LegendrePoly.h"
#include "core/ObjectCache.h"
#include "utils/Printer.h"

using namespace std;
using namespace Eigen;

namespace mrcpp {

typedef ObjectCache<LegendrePoly> LegendreCache;

/** Legendre polynomial constructed on [-1,1] and
  * scaled by n and translated by l */
LegendrePoly::LegendrePoly(int k, double n, double l) :
    Polynomial(k) {
    // Since we create Legendre polynomials recursively on [-1,1]
    // we cache all lower order polynomilas for future use.
    LegendreCache &Cache = LegendreCache::getInstance();
    if (k >= 1) {
        if (not Cache.hasId(k - 1)) {
            LegendrePoly *lp = new LegendrePoly(k - 1);
            Cache.load(k - 1, lp, 2 * sizeof(double) * (k + 1));
        }
    }
    computeLegendrePolynomial(k);
    double a = -1.0;
    double b = 1.0;
    setBounds(&a, &b);
    rescale(n, l);
}

/** Compute Legendre polynomial coefs on interval [-1,1] */
void LegendrePoly::computeLegendrePolynomial(int k) {
    assert(this->size() >= k);
    if (k == 0) {
        this->coefs[0] = 1.0;
    } else if (k == 1) {
        this->coefs[0] = 0.0;
        this->coefs[1] = 1.0;
    } else {
        LegendreCache &Cache = LegendreCache::getInstance();
        LegendrePoly &Lm1 = Cache.get(k - 1);
        LegendrePoly &Lm2 = Cache.get(k - 2);

        double K = (double) k;
        double cm2_0 = Lm2.getCoefs()[0];
        this->coefs[0] = -(K - 1.0)*cm2_0/K;
        for (int j = 1; j < k + 1; j++) {
            double cm1_jm1 = Lm1.getCoefs()[j-1];
            if (j <= k - 2) {
                double cm2_j = Lm2.getCoefs()[j];
                this->coefs[j] = (2.0*K - 1.0)*cm1_jm1/K - (K - 1.0)*cm2_j/K;
            } else {
                this->coefs[j] = (2.0*K - 1.0)*cm1_jm1/K;
            }
        }
    }
}

/** Calculate the value of an n:th order Legendre polynominal in x, including
 * the first derivative.
 */
Vector2d LegendrePoly::firstDerivative(double x) const {
    double c1, c2, c4, ym, yp, y;
    double dy, dyp, dym;

    if (outOfBounds(&x)) {
        MSG_FATAL("Argument out of bounds: " << x << " [" <<
                  this->A[0] << ", " << this->B[0] << "]");
    }

    double q = this->N * x + this->L;
    Vector2d val;

    int order = getOrder();
    if (order == 0) {
        val(0) = 1.0;
        val(1) = 0.0;
        return val;
    }

    if (order == 1) {
        val(0) = q;
        val(1) = this->N * 1.0 + this->L;
        return val;
    }

    y = q;
    dy = 1.0;
    yp = 1.0;
    dyp = 0.0;
    for (int i = 2; i < order + 1; i++) {
        c1 = (double) i;
        c2 = c1 * 2.0 - 1.0;
        c4 = c1 - 1.0;
        ym = y;
        y = (c2 * q * y - c4 * yp) / c1;
        yp = ym;
        dym = dy;
        dy = (c2 * q * dy - c4 * dyp + c2 * yp) / c1;
        dyp = dym;
    }

    val(0) = y;
    val(1) = dy;
    return val;
}

/** Calculate the value of an n:th order Legendre polynominal in x, including
 * first and second derivatives.
 */
Vector3d LegendrePoly::secondDerivative(double x) const {
    NOT_IMPLEMENTED_ABORT;
    double c1, c2, c4, ym, yp, y, d2y;
    double dy, dyp, dym, d2ym, d2yp;

    double q = this->N * x + this->L;
    if (outOfBounds(&x)) {
        MSG_FATAL("Argument out of bounds: " << x << " [" <<
                  this->A[0] << ", " << this->B[0] << "]");
    }

    Vector3d val;

    int order = getOrder();
    if (order == 0) {
        val(0) = 1.e0;
        val(1) = 0.e0;
        val(2) = 0.e0;
        return val;
    }

    if (order == 1) {
        val(0) = q;
        val(1) = this->N * 1.e0 + this->L;
        val(2) = 0.e0;
        return val;
    }

    y = q;
    dy = 1.e0;
    d2y = 0.e0;
    yp = 1.e0;
    dyp = 0.e0;
    d2yp = 0.e0;
    for (int i = 2; i < order + 1; i++) {
        c1 = (double) i;
        c2 = c1 * 2.e0 - 1.e0;
        c4 = c1 - 1.e0;
        ym = y;
        y = (c2 * x * y - c4 * yp) / c1;
        yp = ym;
        dym = dy;
        dy = (c2 * x * dy - c4 * dyp + c2 * yp) / c1;
        dyp = dym;
        d2ym = d2y;
        d2y = (c2 * x * d2y - c4 * d2yp + c2 * 2.e0 * dyp) / c1;
        d2yp = d2ym;
    }
    val(0) = y;
    val(1) = dy;
    val(2) = d2y;
    return val;
}

} // namespace mrcpp
