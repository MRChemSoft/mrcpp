/*
 *
 *  Base class of functions that is representable in the mw basis.
 * This includes gaussians, expansions, polynomials and even function trees.
 */

#pragma once

#include <iostream>
#include <array>

#include "constants.h"

namespace mrcpp {

template<int D>
class RepresentableFunction {
public:
    RepresentableFunction(const double *a = nullptr, const double *b = nullptr);
    RepresentableFunction(const RepresentableFunction<D> &func);
    RepresentableFunction<D> &operator=(const RepresentableFunction<D> &func);
    virtual ~RepresentableFunction();

    virtual double evalf(const double *r) const = 0;
    virtual double evalf(const std::array<double, D> &r) const = 0;

    void setBounds(const double *a, const double *b);
    void clearBounds();

    bool isBounded() const { return this->bounded; }
    bool outOfBounds(const double *r) const;

    double getLowerBound(int d) const { return this->A[d]; }
    double getUpperBound(int d) const { return this->B[d]; }

    const double *getLowerBounds() const { return this->A; }
    const double *getUpperBounds() const { return this->B; }

    virtual bool isVisibleAtScale(int scale, int nQuadPts) const { return true; }
    virtual bool isZeroOnInterval(const double *a, const double *b) const { return false; }

    friend std::ostream& operator<<(std::ostream &o, const RepresentableFunction<D> &func) { return func.print(o); }

protected:
    bool bounded;
    double *A; ///< Lower bound, NULL if unbounded
    double *B; ///< Upper bound, Null if unbounded

    std::ostream& print(std::ostream &o) const;
};

} // namespace mrcpp
