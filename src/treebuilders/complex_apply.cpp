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
#include "apply.h"
#include "ConvolutionCalculator.h"
#include "CopyAdaptor.h"
#include "DefaultCalculator.h"
#include "DerivativeCalculator.h"
#include "SplitAdaptor.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "grid.h"
#include "operators/ConvolutionOperator.h"
#include "operators/DerivativeOperator.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {


/** @brief Application of MW integral convolution operator (complex version)
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] oper: Convolution operator to apply
 * @param[in] inp: Input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree, default -1
 * @param[in] absPrec: Build output tree based on absolute precision, default false
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 * 
 * The default is to work with relative precision
 * (stop when the wavelet coefficients are below a given (small) fraction of
 * function norm.
 * Sometimes it is better to use absolute precision (e.g. a contribution in a sum)
 * which means stop once wavelet coefficients are below a certain (absoute) value
 * Rel prec ∣d∣<ϵ/∣f∣
 * Abs prec ∣d∣<ϵ
 * The two ϵ are not necessarily the same.
 * The first one is (in general) the overall precision of the calculation (not always...)
 * The second one depends on the particular operation which you are performing.
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 * \todo !!! Here should be given a method for greed cleaning !!!
 * 
 *
 */
template <int D>
void apply
(
    double prec, ComplexObject< FunctionTree<D> > &out,
    ComplexObject< ConvolutionOperator<D> > &oper, ComplexObject< FunctionTree<D> > &inp,
    int maxIter, bool absPrec
)
{
    FunctionTree<D> temp1( inp.real->getMRA() );
    FunctionTree<D> temp2( inp.real->getMRA() );

    apply(prec, temp1, *oper.real, *inp.real, maxIter, absPrec);
    apply(prec, temp2, *oper.imaginary, *inp.imaginary, maxIter, absPrec);
    add(prec, *out.real, 1.0, temp1, -1.0, temp2);

    //temp1.setZero();
    //temp2.setZero();

    apply(prec, temp1, *oper.imaginary, *inp.real, maxIter, absPrec);
    apply(prec, temp2, *oper.real, *inp.imaginary, maxIter, absPrec);
    add(prec, *out.imaginary, 1.0, temp1, 1.0, temp2);
}





template
void apply <1>
(
    double prec, ComplexObject< FunctionTree<1> > &out,
    ComplexObject< ConvolutionOperator<1> > &oper, ComplexObject< FunctionTree<1> > &inp,
    int maxIter, bool absPrec
);


} // namespace mrcpp
