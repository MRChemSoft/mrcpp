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
        : coef(c)
        , alpha(a)
        , pos(r) {
}

template<> double Slater<1>::calcSquareNorm() const {
    double c2 = this->coef * this->coef;
    return c2 / this->alpha;
}

template<> double Slater<2>::calcSquareNorm() const {
    double c2 = this->coef * this->coef;
    double k2 = this->alpha * this->alpha;
    return 0.5 * pi * c2 / k2;
}
    
template<> double Slater<3>::calcSquareNorm() const {
    double c2 = this->coef * this->coef;
    double k3 = this->alpha * this->alpha * this->alpha;
    return pi * c2 / k3;
}
    
template <int D> double Slater<D>::evalf(const Coord<D> &r) const {
    auto dist = math_utils::calc_distance<D>(r , pos);
    double exp_val = std::exp(-this->alpha * dist);
    return this->coef * exp_val;
}

template class Slater<1>;
template class Slater<2>;
template class Slater<3>;

} // namespace mrcpp
