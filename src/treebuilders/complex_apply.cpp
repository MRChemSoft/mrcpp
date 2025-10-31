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

#include "complex_apply.h"
#include "ConvolutionCalculator.h"
#include "CopyAdaptor.h"
#include "DefaultCalculator.h"
#include "DerivativeCalculator.h"
#include "SplitAdaptor.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "apply.h"
#include "grid.h"
#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

template <int D>
void apply(double prec,
           ComplexObject<FunctionTree<D>> &out,
           ComplexObject<ConvolutionOperator<D>> &oper,
           ComplexObject<FunctionTree<D>> &inp,
           int maxIter,
           bool absPrec) {
    FunctionTree<D> temp1(inp.real->getMRA());
    FunctionTree<D> temp2(inp.real->getMRA());

    // Real part:  OR*FR - OI*FI
    apply(prec, temp1, *oper.real,      *inp.real,      maxIter, absPrec);
    apply(prec, temp2, *oper.imaginary, *inp.imaginary, maxIter, absPrec);
    add(prec, *out.real, 1.0, temp1, -1.0, temp2);

    // Imag part:  OI*FR + OR*FI
    apply(prec, temp1, *oper.imaginary, *inp.real,      maxIter, absPrec);
    apply(prec, temp2, *oper.real,      *inp.imaginary, maxIter, absPrec);
    add(prec, *out.imaginary, 1.0, temp1, 1.0, temp2);
}

template void apply<1>(double prec,
                       ComplexObject<FunctionTree<1>> &out,
                       ComplexObject<ConvolutionOperator<1>> &oper,
                       ComplexObject<FunctionTree<1>> &inp,
                       int maxIter,
                       bool absPrec);

} // namespace mrcpp