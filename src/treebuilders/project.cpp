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

#include "project.h"
#include "MRCPP/functions/AnalyticFunction.h"
#include "MRCPP/trees/FunctionTree.h"
#include "MRCPP/trees/MultiResolutionAnalysis.h"
#include "MRCPP/utils/Printer.h"
#include "MRCPP/utils/Timer.h"
#include "ProjectionCalculator.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"

namespace mrcpp {

/** @brief Projection of analytic function into MW representation
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp Input function
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed using the general algorithm:
 *  1) Compute MW coefs on current grid
 *  2) Refine grid where necessary based on prec
 *  3) Repeat until convergence or maxIter is reached
 *
 * This algorithm will start at whatever grid is present in the output tree when
 * the function is called (this grid should however be EMPTY, e.i. no coefs).
 *
 * A negative precision means NO refinement, as do maxIter = 0.
 * A negative maxIter means no bound.
 *
 */
template <int D>
void project(double prec, FunctionTree<D> &out, std::function<double(const Coord<D> &r)> func, int maxIter) {
    AnalyticFunction<D> inp(func);
    mrcpp::project(prec, out, inp, maxIter);
}

/** @brief Projection of analytic function into MW representation
 *
 * @param[in] prec Build precision of output function
 * @param[in,out] out Output function to be built
 * @param[in] inp Input function
 * @param[in] maxIter Maximum number of refinement iterations in output tree
 *
 * The output function will be computed using the general algorithm:
 *  1) Compute MW coefs on current grid
 *  2) Refine grid where necessary based on prec
 *  3) Repeat until convergence or maxIter is reached
 *
 * This algorithm will start at whatever grid is present in the output tree when
 * the function is called (this grid should however be EMPTY, e.i. no coefs).
 *
 * A negative precision means NO refinement, as do maxIter = 0.
 * A negative maxIter means no bound.
 *
 */
template <int D> void project(double prec, FunctionTree<D> &out, RepresentableFunction<D> &inp, int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    const auto scaling_factor = out.getMRA().getWorldBox().getScalingFactor();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    ProjectionCalculator<D> calculator(inp, scaling_factor);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Printer::printTime(10, "Time transform", trans_t);
    Printer::printSeparator(10, ' ');
}

template void project<1>(double prec, FunctionTree<1> &out, RepresentableFunction<1> &inp, int maxIter);
template void project<2>(double prec, FunctionTree<2> &out, RepresentableFunction<2> &inp, int maxIter);
template void project<3>(double prec, FunctionTree<3> &out, RepresentableFunction<3> &inp, int maxIter);
template void project<1>(double prec, FunctionTree<1> &out, std::function<double(const Coord<1> &r)> func, int maxIter);
template void project<2>(double prec, FunctionTree<2> &out, std::function<double(const Coord<2> &r)> func, int maxIter);
template void project<3>(double prec, FunctionTree<3> &out, std::function<double(const Coord<3> &r)> func, int maxIter);

} // namespace mrcpp
