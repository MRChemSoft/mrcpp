/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

/*
 *
 *  Base class of functions that is representable in the mw basis.
 * This includes gaussians, expansions, polynomials and even function trees.
 */

#pragma once

#include <iostream>

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

template <int D> class RepresentableFunction {
public:
    RepresentableFunction(const double *a = nullptr, const double *b = nullptr);
    RepresentableFunction(const RepresentableFunction<D> &func);
    RepresentableFunction<D> &operator=(const RepresentableFunction<D> &func);
    virtual ~RepresentableFunction();

    virtual double evalf(const Coord<D> &r) const = 0;

    void setBounds(const double *a, const double *b);
    void clearBounds();

    bool isBounded() const { return this->bounded; }
    bool outOfBounds(const Coord<D> &r) const;

    double getLowerBound(int d) const { return this->A[d]; }
    double getUpperBound(int d) const { return this->B[d]; }

    const double *getLowerBounds() const { return this->A; }
    const double *getUpperBounds() const { return this->B; }

    virtual bool isVisibleAtScale(int scale, int nQuadPts) const { return true; }
    virtual bool isZeroOnInterval(const double *a, const double *b) const { return false; }

    friend std::ostream &operator<<(std::ostream &o, const RepresentableFunction<D> &func) { return func.print(o); }

protected:
    bool bounded;
    double *A; ///< Lower bound, NULL if unbounded
    double *B; ///< Upper bound, Null if unbounded

    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp
