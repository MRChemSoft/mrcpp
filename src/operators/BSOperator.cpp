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

#include "BSOperator.h"
#include "treebuilders/BSCalculator.h"
#include "treebuilders/BandWidthAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @returns New BSOperator object
 *  @param[in] mra: Which MRA the operator is defined
 *  @param[in] order: Derivative order, defined for 1, 2 and 3
 */
template <int D>
BSOperator<D>::BSOperator(const MultiResolutionAnalysis<D> &mra, int order)
        : DerivativeOperator<D>(mra, mra.getRootScale()) {
    this->order = order;
    initialize();
}

template <int D> void BSOperator<D>::initialize() {
    int bw = 1; // Operator bandwidth
    auto oper_mra = this->getOperatorMRA();

    TreeBuilder<2> builder;
    BSCalculator calculator(oper_mra.getScalingBasis(), this->order);
    BandWidthAdaptor adaptor(bw, oper_mra.getMaxScale());

    auto o_tree = std::make_unique<OperatorTree>(oper_mra, MachineZero);
    builder.build(*o_tree, calculator, adaptor, -1);

    Timer trans_t;
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    this->init(1);
    for (int d = 0; d < D; d++) this->assign(0, d, this->raw_exp[0].get());
}

template class BSOperator<1>;
template class BSOperator<2>;
template class BSOperator<3>;

} // namespace mrcpp
