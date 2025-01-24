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

#include "project.h"
#include "ProjectionCalculator.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "functions/AnalyticFunction.h"
#include "trees/FunctionTree.h"
#include "trees/MultiResolutionAnalysis.h"

#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @brief Project an analytic function onto the MW basis, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] inp: Input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D, typename T> void project(double prec, FunctionTree<D, T> &out, std::function<T(const Coord<D> &r)> func, int maxIter, bool absPrec) {
    AnalyticFunction<D, T> inp(func);

    mrcpp::project(prec, out, inp, maxIter, absPrec);
}

/** @brief Project an analytic function onto the MW basis, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] inp: Input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D, typename T> void project(double prec, FunctionTree<D, T> &out, RepresentableFunction<D, T> &inp, int maxIter, bool absPrec) {
    int maxScale = out.getMRA().getMaxScale();
    const auto scaling_factor = out.getMRA().getWorldBox().getScalingFactors();
    TreeBuilder<D, T> builder;
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);

    ProjectionCalculator<D, T> calculator(inp, scaling_factor);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');
}

/** @brief Project an analytic vector function onto the MW basis, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function vector to be built
 * @param[in] inp: Input function vector
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 *
 * @details The output function will be computed using the general algorithm:
 * - Compute MW coefs on current grid
 * - Refine grid where necessary based on `prec`
 * - Repeat until convergence or `maxIter` is reached
 * - `prec < 0` or `maxIter = 0` means NO refinement
 * - `maxIter < 0` means no bound
 *
 * @note This algorithm will start at whatever grid is present in the `out`
 * tree when the function is called (this grid should however be EMPTY, e.i.
 * no coefs).
 *
 */
template <int D, typename T> void project(double prec, FunctionTreeVector<D, T> &out, std::vector<std::function<T(const Coord<D> &r)>> func, int maxIter, bool absPrec) {
    if (out.size() != func.size()) MSG_ABORT("Size mismatch");
    for (auto j = 0; j < D; j++) mrcpp::project<D>(prec, get_func(out, j), func[j], maxIter, absPrec);
}

template void project<1, double>(double prec, FunctionTree<1, double> &out, RepresentableFunction<1, double> &inp, int maxIter, bool absPrec);
template void project<2, double>(double prec, FunctionTree<2, double> &out, RepresentableFunction<2, double> &inp, int maxIter, bool absPrec);
template void project<3, double>(double prec, FunctionTree<3, double> &out, RepresentableFunction<3, double> &inp, int maxIter, bool absPrec);

template void project<1, double>(double prec, FunctionTree<1, double> &out, std::function<double(const Coord<1> &r)> func, int maxIter, bool absPrec);
template void project<2, double>(double prec, FunctionTree<2, double> &out, std::function<double(const Coord<2> &r)> func, int maxIter, bool absPrec);
template void project<3, double>(double prec, FunctionTree<3, double> &out, std::function<double(const Coord<3> &r)> func, int maxIter, bool absPrec);
template void project<1, double>(double prec, FunctionTreeVector<1, double> &out, std::vector<std::function<double(const Coord<1> &r)>> inp, int maxIter, bool absPrec);
template void project<2, double>(double prec, FunctionTreeVector<2, double> &out, std::vector<std::function<double(const Coord<2> &r)>> inp, int maxIter, bool absPrec);
template void project<3, double>(double prec, FunctionTreeVector<3, double> &out, std::vector<std::function<double(const Coord<3> &r)>> inp, int maxIter, bool absPrec);

template void project<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, RepresentableFunction<1, ComplexDouble> &inp, int maxIter, bool absPrec);
template void project<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, RepresentableFunction<2, ComplexDouble> &inp, int maxIter, bool absPrec);
template void project<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, RepresentableFunction<3, ComplexDouble> &inp, int maxIter, bool absPrec);

template void project<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, std::function<ComplexDouble(const Coord<1> &r)> func, int maxIter, bool absPrec);
template void project<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, std::function<ComplexDouble(const Coord<2> &r)> func, int maxIter, bool absPrec);
template void project<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, std::function<ComplexDouble(const Coord<3> &r)> func, int maxIter, bool absPrec);
template void project<1, ComplexDouble>(double prec, FunctionTreeVector<1, ComplexDouble> &out, std::vector<std::function<ComplexDouble(const Coord<1> &r)>> inp, int maxIter, bool absPrec);
template void project<2, ComplexDouble>(double prec, FunctionTreeVector<2, ComplexDouble> &out, std::vector<std::function<ComplexDouble(const Coord<2> &r)>> inp, int maxIter, bool absPrec);
template void project<3, ComplexDouble>(double prec, FunctionTreeVector<3, ComplexDouble> &out, std::vector<std::function<ComplexDouble(const Coord<3> &r)>> inp, int maxIter, bool absPrec);

} // namespace mrcpp
