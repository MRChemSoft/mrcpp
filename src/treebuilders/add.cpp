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

#include <tuple>
#include <vector>

#include "AdditionCalculator.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @brief Addition of two MW function representations, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] a: Numerical coefficient of function a
 * @param[in] inp_a: Input function a
 * @param[in] b: Numerical coefficient of function b
 * @param[in] inp_b: Input function b
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 *
 * @details The output function will be computed as the sum of the two input
 * functions (including the numerical coefficient), using the general algorithm:
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
template <int D>
void add(double prec,
         FunctionTree<D> &out,
         double a,
         FunctionTree<D> &inp_a,
         double b,
         FunctionTree<D> &inp_b,
         int maxIter,
         bool absPrec) {
    FunctionTreeVector<D> tmp_vec;
    tmp_vec.push_back(std::make_tuple(a, &inp_a));
    tmp_vec.push_back(std::make_tuple(b, &inp_b));
    add(prec, out, tmp_vec, maxIter, absPrec);
}

/** @brief Addition of several MW function representations, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] inp: Vector of input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 *
 * @details The output function will be computed as the sum of all input
 * functions in the vector (including their numerical coefficients), using
 * the general algorithm:
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
template <int D> void add(double prec, FunctionTree<D> &out, FunctionTreeVector<D> &inp, int maxIter, bool absPrec) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale, absPrec);
    AdditionCalculator<D> calculator(inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    for (int i = 0; i < inp.size(); i++) {
        FunctionTree<D> &tree = get_func(inp, i);
        tree.deleteGenerated();
    }
    clean_t.stop();

    print::time(10, "Time transform", trans_t);
    print::time(10, "Time cleaning", clean_t);
    print::separator(10, ' ');
}

template void add(double prec,
                  FunctionTree<1> &out,
                  double a,
                  FunctionTree<1> &tree_a,
                  double b,
                  FunctionTree<1> &tree_b,
                  int maxIter,
                  bool absPrec);
template void add(double prec,
                  FunctionTree<2> &out,
                  double a,
                  FunctionTree<2> &tree_a,
                  double b,
                  FunctionTree<2> &tree_b,
                  int maxIter,
                  bool absPrec);
template void add(double prec,
                  FunctionTree<3> &out,
                  double a,
                  FunctionTree<3> &tree_a,
                  double b,
                  FunctionTree<3> &tree_b,
                  int maxIter,
                  bool absPrec);

template void add(double prec, FunctionTree<1> &out, FunctionTreeVector<1> &inp, int maxIter, bool absPrec);
template void add(double prec, FunctionTree<2> &out, FunctionTreeVector<2> &inp, int maxIter, bool absPrec);
template void add(double prec, FunctionTree<3> &out, FunctionTreeVector<3> &inp, int maxIter, bool absPrec);

} // namespace mrcpp
