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

#include <Eigen/Core>

#include "MultiplicationAdaptor.h"
#include "MultiplicationCalculator.h"
#include "PowerCalculator.h"
#include "SquareCalculator.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "add.h"
#include "grid.h"
#include "multiply.h"

#include "trees/FunctionNode.h"
#include "trees/FunctionTree.h"
#include "trees/HilbertIterator.h"

#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/** @brief Multiplication of two MW function representations, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] c: Numerical coefficient
 * @param[in] inp_a: Input function a
 * @param[in] inp_b: Input function b
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 * @param[in] useMaxNorms: Build output tree based on norm estimates from input
 *
 * @details The output function will be computed as the product of the two input
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
void multiply(double prec,
              FunctionTree<D> &out,
              double c,
              FunctionTree<D> &inp_a,
              FunctionTree<D> &inp_b,
              int maxIter,
              bool absPrec,
              bool useMaxNorms) {
    FunctionTreeVector<D> tmp_vec;
    tmp_vec.push_back(std::make_tuple(c, &inp_a));
    tmp_vec.push_back(std::make_tuple(1.0, &inp_b));
    multiply(prec, out, tmp_vec, maxIter, absPrec, useMaxNorms);
}

/** @brief Multiplication of several MW function representations, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] inp: Vector of input function
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 * @param[in] useMaxNorms: Build output tree based on norm estimates from input
 *
 * @details The output function will be computed as the product of all input
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
template <int D>
void multiply(double prec,
              FunctionTree<D> &out,
              FunctionTreeVector<D> &inp,
              int maxIter,
              bool absPrec,
              bool useMaxNorms) {
    for (auto i = 0; i < inp.size(); i++)
        if (out.getMRA() != get_func(inp, i).getMRA()) MSG_ABORT("Incompatible MRA");

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    MultiplicationCalculator<D> calculator(inp);

    if (useMaxNorms) {
        for (int i = 0; i < inp.size(); i++) get_func(inp, i).makeMaxSquareNorms();
        MultiplicationAdaptor<D> adaptor(prec, maxScale, inp);
        builder.build(out, calculator, adaptor, maxIter);
    } else {
        WaveletAdaptor<D> adaptor(prec, maxScale, absPrec);
        builder.build(out, calculator, adaptor, maxIter);
    }

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

/** @brief Out-of-place square of MW function representations, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] inp: Input function to square
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 *
 * @details The output function will be computed as the square of the input
 * function, using the general algorithm:
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
template <int D> void square(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, int maxIter, bool absPrec) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale, absPrec);
    SquareCalculator<D> calculator(inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    inp.deleteGenerated();
    clean_t.stop();

    print::time(10, "Time transform", trans_t);
    print::time(10, "Time cleaning", clean_t);
    print::separator(10, ' ');
}

/** @brief Out-of-place power of MW function representations, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] inp: Input function to square
 * @param[in] p: Numerical power
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 *
 * @details The output function will be computed as the input function raised
 * to the given power, using the general algorithm:
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
void power(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, double p, int maxIter, bool absPrec) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale, absPrec);
    PowerCalculator<D> calculator(inp, p);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    inp.deleteGenerated();
    clean_t.stop();

    print::time(10, "Time transform", trans_t);
    print::time(10, "Time cleaning", clean_t);
    print::separator(10, ' ');
}

/** @brief Dot product of two MW function vectors, adaptive grid
 *
 * @param[in] prec: Build precision of output function
 * @param[out] out: Output function to be built
 * @param[in] inp_a: Input function vector
 * @param[in] inp_b: Input function vector
 * @param[in] maxIter: Maximum number of refinement iterations in output tree
 * @param[in] absPrec: Build output tree based on absolute precision
 *
 * @details The output function will be computed as the dot product of the two
 * input vectors (including their numerical coefficients). The precision
 * parameter is used only in the multiplication part, the final addition will
 * be on the fixed union grid of the components.
 *
 * @note The length of the input vectors must be the same.
 *
 */
template <int D>
void dot(double prec,
         FunctionTree<D> &out,
         FunctionTreeVector<D> &inp_a,
         FunctionTreeVector<D> &inp_b,
         int maxIter,
         bool absPrec) {
    if (inp_a.size() != inp_b.size()) MSG_ABORT("Input length mismatch");

    FunctionTreeVector<D> tmp_vec;
    for (int d = 0; d < inp_a.size(); d++) {
        double coef_a = get_coef(inp_a, d);
        double coef_b = get_coef(inp_b, d);
        FunctionTree<D> &tree_a = get_func(inp_a, d);
        FunctionTree<D> &tree_b = get_func(inp_b, d);
        auto *out_d = new FunctionTree<D>(out.getMRA());
        build_grid(*out_d, out);
        multiply(prec, *out_d, 1.0, tree_a, tree_b, maxIter, absPrec);
        tmp_vec.push_back(std::make_tuple(coef_a * coef_b, out_d));
    }
    build_grid(out, tmp_vec);
    add(-1.0, out, tmp_vec, 0);
    clear(tmp_vec, true);
}

/** @returns Dot product <bra|ket> of two MW function representations
 *
 * @param[in] bra: Bra side input function
 * @param[in] ket: Ket side input function
 *
 * @details The dot product is computed with the trees in compressed form, i.e.
 * scaling coefs only on root nodes, wavelet coefs on all nodes. Since wavelet
 * functions are orthonormal through ALL scales and the root scaling functions
 * are orthonormal to all finer level wavelet functions, this becomes a rather
 * efficient procedure as you only need to compute the dot product where the
 * grids overlap.
 *
 */
template <int D> double dot(FunctionTree<D> &bra, FunctionTree<D> &ket) {
    if (bra.getMRA() != ket.getMRA()) MSG_ABORT("Trees not compatible");

    MWNodeVector<D> nodeTable;
    HilbertIterator<D> it(&bra);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        nodeTable.push_back(&node);
    }
    int nNodes = nodeTable.size();
    double result = 0.0;
    double locResult = 0.0;
    // OMP is disabled in order to get EXACT results (to the very last digit), the
    // order of summation makes the result different beyond the 14th digit or so.
    // OMP does improve the performace, but its not worth it for the time being.
    //#pragma omp parallel firstprivate(n_nodes, locResult) num_threads(mrcpp_get_num_threads())
    //		shared(nodeTable,rhs,result)
    //    {
    //#pragma omp for schedule(guided)
    for (int n = 0; n < nNodes; n++) {
        const auto &braNode = static_cast<const FunctionNode<D> &>(*nodeTable[n]);
        const MWNode<D> *mwNode = ket.findNode(braNode.getNodeIndex());
        if (mwNode == nullptr) continue;

        const auto &ketNode = static_cast<const FunctionNode<D> &>(*mwNode);
        if (braNode.isRootNode()) { locResult += dotScaling(braNode, ketNode); }
        locResult += dotWavelet(braNode, ketNode);
    }
    //#pragma omp critical
    result += locResult;
    //    }
    return result;
}

/** @brief abs-dot product of two MW function representations
 *
 * @param[in] bra: Bra side input function
 * @param[in] ket: Ket side input function
 *
 * If exact=true: the grid of ket MUST include the grid of bra.
 * If exact=false: does not at any time read the coefficients individually.
 * The product is done for the end nodes of the bra multiplied by the nodes from the
 * ket with either the same idx, or using a lower scale and assuming uniform
 * distribution within the node.
 * If the product is zero, the functions are disjoints.
 */
template <int D> double node_norm_dot(FunctionTree<D> &bra, FunctionTree<D> &ket, bool exact) {
    if (bra.getMRA() != ket.getMRA()) MSG_ABORT("Incompatible MRA");

    double result = 0.0;
    int ncoef = bra.getKp1_d() * bra.getTDim();
    double valA[ncoef];
    double valB[ncoef];
    int nNodes = bra.getNEndNodes();

    for (int n = 0; n < nNodes; n++) {
        FunctionNode<D> &node = bra.getEndFuncNode(n);
        const NodeIndex<D> idx = node.getNodeIndex();
        if (exact) {
            // convert to interpolating coef, take abs, convert back
            FunctionNode<D> *mwNode = static_cast<FunctionNode<D> *>(ket.findNode(idx));
            if (mwNode == nullptr) MSG_ABORT("Trees must have same grid");
            node.getAbsCoefs(valA);
            mwNode->getAbsCoefs(valB);
            for (int i = 0; i < ncoef; i++) result += valA[i] * valB[i];
        } else {
            // approximate by product of node norms
            int rIdx = ket.getRootBox().getBoxIndex(idx);
            assert(rIdx >= 0);
            const MWNode<D> &root = ket.getRootBox().getNode(rIdx);
            result += std::sqrt(node.getSquareNorm()) * root.getNodeNorm(idx);
        }
    }

    return result;
}

template void multiply<1>(double prec,
                          FunctionTree<1> &out,
                          double c,
                          FunctionTree<1> &tree_a,
                          FunctionTree<1> &tree_b,
                          int maxIter,
                          bool absPrec,
                          bool useMaxNorms);
template void multiply<2>(double prec,
                          FunctionTree<2> &out,
                          double c,
                          FunctionTree<2> &tree_a,
                          FunctionTree<2> &tree_b,
                          int maxIter,
                          bool absPrec,
                          bool useMaxNorms);
template void multiply<3>(double prec,
                          FunctionTree<3> &out,
                          double c,
                          FunctionTree<3> &tree_a,
                          FunctionTree<3> &tree_b,
                          int maxIter,
                          bool absPrec,
                          bool useMaxNorms);
template void multiply<1>(double prec,
                          FunctionTree<1> &out,
                          FunctionTreeVector<1> &inp,
                          int maxIter,
                          bool absPrec,
                          bool useMaxNorms);
template void multiply<2>(double prec,
                          FunctionTree<2> &out,
                          FunctionTreeVector<2> &inp,
                          int maxIter,
                          bool absPrec,
                          bool useMaxNorms);
template void multiply<3>(double prec,
                          FunctionTree<3> &out,
                          FunctionTreeVector<3> &inp,
                          int maxIter,
                          bool absPrec,
                          bool useMaxNorms);
template void power<1>(double prec, FunctionTree<1> &out, FunctionTree<1> &tree, double pow, int maxIter, bool absPrec);
template void power<2>(double prec, FunctionTree<2> &out, FunctionTree<2> &tree, double pow, int maxIter, bool absPrec);
template void power<3>(double prec, FunctionTree<3> &out, FunctionTree<3> &tree, double pow, int maxIter, bool absPrec);
template void square<1>(double prec, FunctionTree<1> &out, FunctionTree<1> &tree, int maxIter, bool absPrec);
template void square<2>(double prec, FunctionTree<2> &out, FunctionTree<2> &tree, int maxIter, bool absPrec);
template void square<3>(double prec, FunctionTree<3> &out, FunctionTree<3> &tree, int maxIter, bool absPrec);
template void dot<1>(double prec,
                     FunctionTree<1> &out,
                     FunctionTreeVector<1> &inp_a,
                     FunctionTreeVector<1> &inp_b,
                     int maxIter,
                     bool absPrec);
template void dot<2>(double prec,
                     FunctionTree<2> &out,
                     FunctionTreeVector<2> &inp_a,
                     FunctionTreeVector<2> &inp_b,
                     int maxIter,
                     bool absPrec);
template void dot<3>(double prec,
                     FunctionTree<3> &out,
                     FunctionTreeVector<3> &inp_a,
                     FunctionTreeVector<3> &inp_b,
                     int maxIter,
                     bool absPrec);
template double dot<1>(FunctionTree<1> &bra, FunctionTree<1> &ket);
template double dot<2>(FunctionTree<2> &bra, FunctionTree<2> &ket);
template double dot<3>(FunctionTree<3> &bra, FunctionTree<3> &ket);

template double node_norm_dot<1>(FunctionTree<1> &bra, FunctionTree<1> &ket, bool exact);
template double node_norm_dot<2>(FunctionTree<2> &bra, FunctionTree<2> &ket, bool exact);
template double node_norm_dot<3>(FunctionTree<3> &bra, FunctionTree<3> &ket, bool exact);

} // namespace mrcpp
