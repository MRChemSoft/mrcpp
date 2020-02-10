/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "ABGVOperator.h"
#include "treebuilders/ABGVCalculator.h"
#include "treebuilders/BandWidthAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "trees/OperatorTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @returns New ABGVOperator object
 *  @param[in] mra: Which MRA the operator is defined
 *  @param[in] a: Left boundary condition
 *  @param[in] b: Right boundary condition
 *  @details Boundary parameters correspond to:
 *  - `a=0.0` `b=0.0`: Strictly local "center" difference
 *  - `a=0.5` `b=0.5`: Semi-local central difference
 *  - `a=1.0` `b=0.0`: Semi-local forward difference
 *  - `a=0.0` `b=1.0`: Semi-local backward difference
 */
template <int D>
ABGVOperator<D>::ABGVOperator(const MultiResolutionAnalysis<D> &mra, double a, double b)
        : DerivativeOperator<D>(mra) {
    initializeOperator(a, b);
}

template <int D> void ABGVOperator<D>::initializeOperator(double a, double b) {
    int bw = 0; // Operator bandwidth
    if (std::abs(a) > MachineZero) bw = 1;
    if (std::abs(b) > MachineZero) bw = 1;
    int max_scale = this->oper_mra.getMaxScale();
    const ScalingBasis &basis = this->oper_mra.getScalingBasis();

    TreeBuilder<2> builder;
    ABGVCalculator calculator(basis, a, b);
    BandWidthAdaptor adaptor(bw, max_scale);

    auto *o_tree = new OperatorTree(this->oper_mra, MachineZero);
    builder.build(*o_tree, calculator, adaptor, -1);

    Timer trans_t;
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->oper_exp.push_back(o_tree);
}

template class ABGVOperator<1>;
template class ABGVOperator<2>;
template class ABGVOperator<3>;

} // namespace mrcpp
