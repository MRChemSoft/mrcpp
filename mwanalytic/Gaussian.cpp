/**
 *
 *
 * \date May 25, 2010
 * \author Stig Rune Jensen
 *		   CTCC, University of Tromsø
 *
 */

#include "Gaussian.h"
#include "NodeIndex.h"

using namespace Eigen;
using namespace std;

template<int D>
Gaussian<D>::Gaussian(double a, double c, const double r[D], const int p[D]) {
    this->alpha = a;
    this->coef = c;
    this->screen = false;
    for (int d = 0; d < D; d++) {
        if (r == 0) {
            this->pos[d] = 0.0;
        } else {
            this->pos[d] = r[d];
        }
        if (p == 0) {
            this->power[d] = 0;
        } else {
            this->power[d] = p[d];
        }
    }
    this->squareNorm = -1.0;
}

template<int D>
Gaussian<D>::~Gaussian() {
}

template<int D>
void Gaussian<D>::multPureGauss(const Gaussian<D> &lhs,	const Gaussian<D> &rhs) {
    double newPos[D], relPos[D];
    double newAlpha = lhs.alpha + rhs.alpha;
    double newCoef = 1.0;
    double mju = (lhs.alpha * rhs.alpha) / newAlpha;
    for (int d = 0; d < D; d++) {
        newPos[d] = (lhs.alpha*lhs.pos[d] + rhs.alpha*rhs.pos[d])/newAlpha;
        relPos[d] = lhs.pos[d] - rhs.pos[d];
        newCoef *= exp(-mju * pow(relPos[d], 2.0));
    }
    setExp(newAlpha);
    setPos(newPos);
    this->squareNorm = -1.0;
    setCoef(newCoef);
}

template<int D>
void Gaussian<D>::calcScreening(double nStdDev) {
    assert(nStdDev > 0);
    double limit = sqrt(nStdDev/this->alpha);
    if (not this->isBounded()) {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
    }
    for (int d = 0; d < D; d++) {
        this->A[d] = this->pos[d] - limit;
        this->B[d] = this->pos[d] + limit;
    }
    screen = true;
}

template<int D>
bool Gaussian<D>::checkScreen(int n, const int *l) const {
    if (not getScreen()) {
        return false;
    }
    double length = pow(2.0, -n);
    const double *A = this->getLowerBounds();
    const double *B = this->getUpperBounds();
    for (int d = 0; d < D; d++) {
        double a = length * l[d];
        double b = length * (l[d] + 1);
        if (a > B[d] or b < A[d]) {
            return true;
        }
    }
    return false;
}

template<int D>
bool Gaussian<D>::isVisibleAtScale(int scale, int nQuadPts) const {
    double stdDeviation = pow(2.0*this->alpha, -0.5);
    int visibleScale = int(-floor(log2(nQuadPts*2.0*stdDeviation)));
    if (scale < visibleScale) {
        return false;
    }
    return true;
}

template<int D>
bool Gaussian<D>::isZeroOnInterval(const double *a, const double *b) const {
    double stdDeviation = pow(2.0*this->alpha, -0.5);
    double gaussBoxMin;
    double gaussBoxMax;
    for (int i=0; i < D; i++) {
        gaussBoxMin = this->pos[i] - 5.0*stdDeviation;
        gaussBoxMax = this->pos[i] + 5.0*stdDeviation;
        if (a[i] > gaussBoxMax or b[i] < gaussBoxMin) {
            return true;
        }
    }
    return false;
}

template<int D>
void Gaussian<D>::evalf(const MatrixXd &points, MatrixXd &values) const {
    assert(points.cols() == D);
    assert(points.cols() == values.cols());
    assert(points.rows() == values.rows());
    for (int d = 0; d < D; d++) {
        for (int i = 0; i < points.rows(); i++) {
            values(i, d) = evalf(points(i, d), d);
        }
    }
}

template class Gaussian<1>;
template class Gaussian<2>;
template class Gaussian<3>;
