/**
 *
 *  Base class for Gaussian type functions
 */

#pragma once

#include <iostream>
#include <cmath>

#pragma GCC system_header
#include <Eigen/Core>

#include "RepresentableFunction.h"
#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class Gaussian: public RepresentableFunction<D> {
public:
    Gaussian(double a, double c, const double r[D], const int p[D]);
    Gaussian(double a, double c, const Coord<D> &r, const std::array<int, D> &p);
    Gaussian<D> &operator=(const Gaussian<D> &gp) = delete;
    virtual Gaussian<D> *copy() const = 0;
    virtual ~Gaussian();

    virtual double evalf(const Coord<D> &r) const = 0;
    virtual double evalf(double r, int dim) const = 0;
    void evalf(const Eigen::MatrixXd &points, Eigen::MatrixXd &values) const;

    virtual double calcSquareNorm() = 0;
    virtual double calcOverlap(GaussFunc<D> &b) = 0;
    virtual double calcOverlap(GaussPoly<D> &b) = 0;

    virtual GaussPoly<D> differentiate(int dir) = 0;

    void calcScreening(double stdDeviations);

    double getSquareNorm() {
        if (this->squareNorm < 0.0) {
            calcSquareNorm();
        }
        return this->squareNorm;
    }
    void normalize() {
        double norm = std::sqrt(getSquareNorm());
        multConstInPlace(1.0/norm);
    }
    void multPureGauss(const Gaussian<D> &lhs, const Gaussian<D> &rhs);
    void multConstInPlace(double c) {
        this->coef *= c;
        this->calcSquareNorm();
    }
    void operator*=(double c) {
        multConstInPlace(c);
    }

    bool checkScreen(int n, const int *l) const;
    int getPower(int i) const { return power[i]; }
    bool getScreen() const { return screen; }
    const int *getPower() const { return power; }
    const double *getPos() const { return pos; }
    double getCoef() const { return coef; }
    auto getExp() const { return alpha[0]; }

    virtual void setPower(const int power[D]) = 0;
    virtual void setPower(int d, int power) = 0;
    void setScreen(bool _screen) { this->screen = _screen; }
    void setCoef(double cf) { this->coef = cf; this->squareNorm = -1.0; }
    void setExp(double _alpha) { this->alpha.fill(_alpha); this->squareNorm = -1.0; }
    void setPos(const double r[D]) {
        for (int d = 0; d < D; d++) {
            this->pos[d] = r[d];
        }
        this->squareNorm = -1.0;
    }
    friend std::ostream& operator<<(std::ostream &o, const Gaussian<D> &gauss) { return gauss.print(o);  }

    friend class GaussExp<D>;
protected:
    bool screen;
    int power[D];		/**< max power in each dim  */
    double coef;		/**< constant factor */
     // double alpha;		/**< exponent  */
    std::array<double, D> alpha;		/**< exponent  */
    double pos[D];		/**< center  */
    double squareNorm;

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

    virtual std::ostream& print(std::ostream &o) const = 0;
};

}
