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

#include "CartesianConvolution.h"

#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"

#include "functions/GaussExp.h"
#include "functions/Gaussian.h"

#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"

#include "trees/BandWidth.h"
#include "trees/FunctionTreeVector.h"
#include "trees/OperatorTree.h"

#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

namespace mrcpp {

CartesianConvolution::CartesianConvolution(const MultiResolutionAnalysis<3> &mra, GaussExp<1> &kernel, double prec)
        : ConvolutionOperator<3>(mra)
        , sep_rank(kernel.size()) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    auto o_prec = prec;
    auto k_prec = prec / 10.0;

    for (auto &k : kernel) k->setPow({0});
    this->initialize(kernel, k_prec, o_prec);
    for (auto &k : kernel) k->setPow({1});
    this->initialize(kernel, k_prec, o_prec);
    for (auto &k : kernel) k->setPow({2});
    this->initialize(kernel, k_prec, o_prec);

    this->initOperExp(this->sep_rank);
    Printer::setPrintLevel(oldlevel);
}

void CartesianConvolution::setCartesianComponents(int x, int y, int z) {
    int x_shift = x * this->sep_rank;
    int y_shift = y * this->sep_rank;
    int z_shift = z * this->sep_rank;

    for (int i = 0; i < this->sep_rank; i++) this->assign(i, 0, this->raw_exp[x_shift + i].get());
    for (int i = 0; i < this->sep_rank; i++) this->assign(i, 1, this->raw_exp[y_shift + i].get());
    for (int i = 0; i < this->sep_rank; i++) this->assign(i, 2, this->raw_exp[z_shift + i].get());
}

} // namespace mrcpp
