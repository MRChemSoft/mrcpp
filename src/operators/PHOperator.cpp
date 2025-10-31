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

#include "PHOperator.h"
#include "treebuilders/BandWidthAdaptor.h"
#include "treebuilders/PHCalculator.h"
#include "treebuilders/TreeBuilder.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

template <int D>
PHOperator<D>::PHOperator(const MultiResolutionAnalysis<D> &mra, int order)
        : DerivativeOperator<D>(mra, mra.getRootScale(), -10) {
    this->order = order;
    initialize();
}

template <int D> void PHOperator<D>::initialize() {
    auto o_mra = this->getOperatorMRA();

    TreeBuilder<2> builder;

    auto &basis = this->MRA.getScalingBasis();
    PHCalculator calculator(basis, this->order);

    int bw = 1;
    int max_scale = this->MRA.getMaxScale();
    BandWidthAdaptor adaptor(bw, max_scale);

    auto o_tree = std::make_unique<OperatorTree>(o_mra, MachineZero);
    builder.build(*o_tree, calculator, adaptor, -1);

    Timer trans_t;
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));
    this->initOperExp(1);
}

template class PHOperator<1>;
template class PHOperator<2>;
template class PHOperator<3>;

} // namespace mrcpp