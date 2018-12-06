/*
 *
 *
 * \date May 25, 2010
 * \author Stig Rune Jensen
 *		   CTCC, University of Tromsø
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include <vector>

#include "Gaussian.h"
#include "Polynomial.h"
#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class GaussPoly final : public Gaussian<D> {
public:
    GaussPoly(double alpha = 0.0, double coef = 1.0, const double pos[D] = nullptr, const int power[D] = nullptr);
    GaussPoly(double alpha, double coef, const Coord<D> &pos, const std::array<int, D> &power);
    GaussPoly(const std::array<double, D> &alpha, double coef, const Coord<D> &pos, const std::array<int, D> &power);
    GaussPoly(const GaussPoly<D> &gp);
    GaussPoly(const GaussFunc<D> &gf);
    GaussPoly<D> &operator=(const GaussPoly<D> &gp) = delete;
    Gaussian<D> *copy() const;
    ~GaussPoly();

    double calcSquareNorm();

    double evalf(const Coord<D> &r) const;
    double evalf(double r, int dim) const;

    double calcOverlap(GaussFunc<D> &b);
    double calcOverlap(GaussPoly<D> &b);

    GaussPoly differentiate(int dir);

    void multInPlace(const GaussPoly<D> &rhs);
    void operator*=(const GaussPoly<D> &rhs) { multInPlace(rhs); }
    GaussPoly<D> mult(const GaussPoly<D> &rhs);
    GaussPoly<D> mult(double c);
    GaussPoly<D> operator*(const GaussPoly<D> &rhs) { return mult(rhs); }
    GaussPoly<D> operator*(double c) { return mult(c); }

    const Eigen::VectorXd &getPolyCoefs(int i) const { return poly[i]->getCoefs(); }
    Eigen::VectorXd &getPolyCoefs(int i) { return poly[i]->getCoefs(); }
    const Polynomial &getPoly(int i) const { return *poly[i]; }
    Polynomial &getPoly(int i) { return *poly[i]; }

    void setPower(int d, int pow);
    void setPower(const std::array<int, D> &pow);
    void setPower(const int pow[D]);
    void setPoly(int d, Polynomial &poly);

    void fillCoefPowVector(std::vector<double> &coefs, std::vector<int *> &power, int pow[D], int dir) const;

private:
    Polynomial *poly[D];

    std::ostream& print(std::ostream &o) const;
};

}
