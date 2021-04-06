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

#include "BoysFunction.h"
#include "core/InterpolatingBasis.h"
#include "treebuilders/project.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"

namespace mrcpp {

BoysFunction::BoysFunction(int n, double p)
        : RepresentableFunction<1>()
        , order(n)
        , prec(p)
        , MRA(BoundingBox<1>(), InterpolatingBasis(13)) {}

double BoysFunction::evalf(const Coord<1> &r) const {
    int oldlevel = Printer::setPrintLevel(0);

    int n = this->order;
    double x = r[0];
    auto f = [x, n](const Coord<1> &t) -> double {
        double t_2 = t[0] * t[0];
        double xt_2 = x * t_2;
        double t_2n = 1.0;
        if (n > 0) { t_2n = std::pow(t_2, n); }
        return std::exp(-xt_2) * t_2n;
    };

    FunctionTree<1> tree(this->MRA);
    mrcpp::project<1>(this->prec, tree, f);
    double result = tree.integrate();

    Printer::setPrintLevel(oldlevel);
    return result;
}

} // namespace mrcpp
