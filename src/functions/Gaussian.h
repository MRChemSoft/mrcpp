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
    Gaussian(double a, double c, const Coord<D> &r, const std::array<int, D> &p);
    Gaussian(const std::array<double, D> &a, double c, const Coord<D> &r, const std::array<int, D> &p);
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
    const std::array<int, D> &getPower() const { return power; }
    const std::array<double, D> &getPos() const { return pos; }
    double getCoef() const { return coef; }
    auto getExp() const { return alpha; }

    virtual void setPower(const std::array<int, D> &power) = 0;
    virtual void setPower(int d, int power) = 0;
    void setScreen(bool _screen) { this->screen = _screen; }
    void setCoef(double cf) { this->coef = cf; this->squareNorm = -1.0; }
    void setExp(double _alpha) { this->alpha.fill(_alpha); this->squareNorm = -1.0; }
    void setExp(const std::array<double, D> &_alpha) { this->alpha = _alpha; this->squareNorm = -1.0; }
    void setPos(const std::array<double, D> &r) { this->pos = r; }

    friend std::ostream& operator<<(std::ostream &o, const Gaussian<D> &gauss) { return gauss.print(o);  }

    friend class GaussExp<D>;
protected:
    bool screen;
    double coef;		            /**< constant factor */
    std::array<int, D> power;		/**< max power in each dim  */
    std::array<double, D> alpha;    /**< exponent  */
    std::array<double, D> pos;		/**< center  */
    double squareNorm;

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

    virtual std::ostream& print(std::ostream &o) const = 0;
};

} // namespace mrcpp
