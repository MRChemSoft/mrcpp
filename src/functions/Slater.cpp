/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include <numeric>

#include "Slater.h"
#include "function_utils.h"
#include "trees/NodeIndex.h"
#include "utils/Printer.h"
#include "utils/details.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

template <int D>
Slater<D>::Slater(double a, double c, const Coord<D> &r)
        : screen(false)
        , coef(c)
        , pos(r) {
        , alpha(a);
}

template <int D> void Slater<D>::calcScreening(double nStdDev) {
    MSG_ABORT("Not Implemented")
}

template <int D> bool Slater<D>::checkScreen(int n, const int *l) const {
    MSG_ABORT("Not Implemented")
}

template <int D> bool Slater<D>::isVisibleAtScale(int scale, int nQuadPts) const {
    MSG_ABORT("Not Implemented")
}

template <int D> bool Slater<D>::isZeroOnInterval(const double *a, const double *b) const {
    MSG_ABORT("Not Implemented")
}

template <int D> void Slater<D>::evalf(const MatrixXd &points, MatrixXd &values) const {
    assert(points.cols() == D);
    assert(points.cols() == values.cols());
    assert(points.rows() == values.rows());
    for (int d = 0; d < D; d++) {
        for (int i = 0; i < points.rows(); i++) { values(i, d) = evalf1D(points(i, d), d); }
    }
}

template class Slater<1>;
template class Slater<2>;
template class Slater<3>;

} // namespace mrcpp
