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

/**
 * @file multiply.cpp
 * @brief Adaptive algebra on multiresolution (MW) function trees: product,
 *        square, power, (componentwise) dot, and inner products.
 *
 * @details
 * This module implements a family of adaptive build routines that produce
 * a new #mrcpp::FunctionTree from algebraic combinations of one or more input
 * trees. The build is driven by the multiresolution refinement loop
 * (TreeBuilder + Adaptor + Calculator):
 *
 *  - On the current output grid, local contributions are computed by a
 *    Calculator (e.g. MultiplicationCalculator, SquareCalculator, PowerCalculator).
 *  - A refinement Adaptor (WaveletAdaptor by default, or MultiplicationAdaptor
 *    when useMaxNorms is enabled) decides whether to split nodes based on
 *    requested precision.
 *  - The refine–recompute process repeats until the target precision is met
 *    or the iteration limit is reached.
 *
 * Precision semantics:
 *  - Relative precision (absPrec = false): split while |d| is not small
 *    relative to the function norm.
 *  - Absolute precision (absPrec = true): split while |d| is above a fixed
 *    absolute threshold.
 *
 * Notes:
 *  - All routines assume the output tree starts with an empty grid (no coeffs).
 *    The grid is grown adaptively unless otherwise stated.
 *  - The input and output trees must belong to compatible MRAs.
 *  - Some routines can optionally use max-norm estimates from inputs to guide
 *    refinement (useMaxNorms).
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
#include "trees/TreeIterator.h"

#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/**
 * @brief Adaptive product of two MW functions with an overall scalar factor.
 *
 * @tparam D  Spatial dimension (1, 2, or 3).
 * @tparam T  Coefficient type (double or ComplexDouble).
 *
 * @param[in]  prec       Target build precision.
 * @param[out] out        Output function tree to construct.
 * @param[in]  c          Scalar prefactor multiplying inp_a * inp_b.
 * @param[in]  inp_a      First input tree.
 * @param[in]  inp_b      Second input tree.
 * @param[in]  maxIter    Max refinement iterations (-1 means unbounded).
 * @param[in]  absPrec    If true: absolute precision; else relative.
 * @param[in]  useMaxNorms If true: use MultiplicationAdaptor with local
 *                         max-norm estimates from inputs for split checks.
 * @param[in]  conjugate  If true: apply complex conjugation to inp_b during multiplication.
 *
 * @details
 * Builds out = c * inp_a * (conjugate ? conj(inp_b) : inp_b) on an adaptively
 * refined grid. If useMaxNorms is true, each input tree contributes local
 * estimates (makeMaxSquareNorms) to scale the precision per node.
 */
template <int D, typename T>
void multiply(double prec, FunctionTree<D, T> &out, T c, FunctionTree<D, T> &inp_a, FunctionTree<D, T> &inp_b, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate) {
    FunctionTreeVector<D, T> tmp_vec;
    tmp_vec.push_back({c, &inp_a});
    tmp_vec.push_back({1.0, &inp_b});
    multiply(prec, out, tmp_vec, maxIter, absPrec, useMaxNorms, conjugate);
}

/**
 * @brief Adaptive product of several MW functions (with per-input scalars).
 *
 * @tparam D  Spatial dimension (1, 2, or 3).
 * @tparam T  Coefficient type.
 *
 * @param[in]  prec        Target build precision.
 * @param[out] out         Output function tree to construct.
 * @param[in]  inp         Vector of inputs (scalar, tree) pairs.
 * @param[in]  maxIter     Max refinement iterations (-1 means unbounded).
 * @param[in]  absPrec     If true: absolute precision; else relative.
 * @param[in]  useMaxNorms Use norm-based adaptor when true.
 * @param[in]  conjugate   Conjugate all trees except the first (if complex).
 *
 * @details
 * Builds out = (Π_k a_k * f_k) where each (a_k, f_k) is the k-th pair.
 * If conjugate is true, all factors except the first are conjugated in the
 * complex case. When useMaxNorms is true, #mrcpp::MultiplicationAdaptor
 * scales the split threshold by input-node max norms to improve targeting.
 */
template <int D, typename T>
void multiply(double prec, FunctionTree<D, T> &out, FunctionTreeVector<D, T> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate) {
    for (auto i = 0; i < inp.size(); i++)
        if (out.getMRA() != get_func(inp, i).getMRA()) MSG_ABORT("Incompatible MRA");

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    MultiplicationCalculator<D, T> calculator(inp, conjugate);

    if (useMaxNorms) {
        for (int i = 0; i < inp.size(); i++) get_func(inp, i).makeMaxSquareNorms();
        MultiplicationAdaptor<D, T> adaptor(prec, maxScale, inp);
        builder.build(out, calculator, adaptor, maxIter);
    } else {
        WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
        builder.build(out, calculator, adaptor, maxIter);
    }

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    for (int i = 0; i < inp.size(); i++) {
        FunctionTree<D, T> &tree = get_func(inp, i);
        tree.deleteGenerated();
    }
    clean_t.stop();

    print::time(10, "Time transform", trans_t);
    print::time(10, "Time cleaning", clean_t);
    print::separator(10, ' ');
}

/**
 * @brief Convenience overload: product of a list of trees (unit coefficients).
 *
 * @tparam D  Spatial dimension.
 * @tparam T  Coefficient type.
 */
template <int D, typename T>
void multiply(double prec, FunctionTree<D, T> &out, std::vector<FunctionTree<D, T> *> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate) {
    FunctionTreeVector<D, T> inp_vec;
    for (auto &t : inp) inp_vec.push_back({1.0, t});
    multiply(prec, out, inp_vec, maxIter, absPrec, useMaxNorms, conjugate);
}

/**
 * @brief Adaptive, out-of-place square: out = (conjugate ? conj(inp) : inp)^2.
 *
 * @tparam D  Spatial dimension.
 * @tparam T  Coefficient type.
 *
 * @details
 * Uses #mrcpp::SquareCalculator over a wavelet-driven adaptive refinement.
 */
template <int D, typename T>
void square(double prec, FunctionTree<D, T> &out, FunctionTree<D, T> &inp, int maxIter, bool absPrec, bool conjugate) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    SquareCalculator<D, T> calculator(inp, conjugate);

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

/**
 * @brief Adaptive power: out = inp^p (real exponent p).
 *
 * @tparam D  Spatial dimension.
 * @tparam T  Coefficient type.
 *
 * @warning Conjugated inputs are not supported here.
 */
template <int D, typename T>
void power(double prec, FunctionTree<D, T> &out, FunctionTree<D, T> &inp, double p, int maxIter, bool absPrec) {
    if (out.getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    if (inp.conjugate()) MSG_ABORT("Not implemented");

    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D, T> builder;
    WaveletAdaptor<D, T> adaptor(prec, maxScale, absPrec);
    PowerCalculator<D, T> calculator(inp, p);

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

/**
 * @brief Adaptive componentwise dot product of two function vectors.
 *
 * @tparam D  Spatial dimension.
 * @tparam T  Coefficient type.
 *
 * @param[in]  prec     Target build precision for the per-component products.
 * @param[out] out      Output tree holding the sum over component products.
 * @param[in]  inp_a    First vector of (scalar, tree) pairs.
 * @param[in]  inp_b    Second vector of (scalar, tree) pairs.
 * @param[in]  maxIter  Max refinement iterations per component product.
 * @param[in]  absPrec  Absolute vs relative precision.
 *
 * @details
 * Computes out = Σ_d (a_d f_d) · (b_d g_d) by first forming per-component
 * products on grids compatible with @p out, then summing these contributions
 * on the fixed union grid (addition step uses a fixed grid, not adaptive).
 */
template <int D, typename T>
void dot(double prec, FunctionTree<D, T> &out, FunctionTreeVector<D, T> &inp_a, FunctionTreeVector<D, T> &inp_b, int maxIter, bool absPrec) {
    if (inp_a.size() != inp_b.size()) MSG_ABORT("Input length mismatch");

    FunctionTreeVector<D, T> tmp_vec;
    for (int d = 0; d < inp_a.size(); d++) {
        T coef_a = get_coef(inp_a, d);
        T coef_b = get_coef(inp_b, d);
        FunctionTree<D, T> &tree_a = get_func(inp_a, d);
        FunctionTree<D, T> &tree_b = get_func(inp_b, d);
        auto *out_d = new FunctionTree<D, T>(out.getMRA());
        build_grid(*out_d, out);
        T One = 1.0;
        multiply(prec, *out_d, One, tree_a, tree_b, maxIter, absPrec, true);
        tmp_vec.push_back({coef_a * coef_b, out_d});
    }
    build_grid(out, tmp_vec);
    add(-1.0, out, tmp_vec, 0);
    clear(tmp_vec, true);
}

/**
 * @brief Inner product ⟨bra|ket⟩ on compressed MW trees.
 *
 * @tparam D  Spatial dimension.
 * @tparam T  Coefficient type of bra.
 * @tparam U  Coefficient type of ket.
 * @tparam V  Return type (double or ComplexDouble).
 *
 * @details
 * Works directly on compressed representation: scaling coefficients on roots
 * and wavelet coefficients on all nodes. Orthonormality across scales makes
 * this efficient: only overlapping nodes contribute.
 */
template <int D, typename T, typename U, typename V>
V dot(FunctionTree<D, T> &bra, FunctionTree<D, U> &ket) {
    if (bra.getMRA() != ket.getMRA()) MSG_ABORT("Trees not compatible");
    MWNodeVector<D, T> nodeTable;
    TreeIterator<D, T> it(bra);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D, T> &node = it.getNode();
        nodeTable.push_back(&node);
    }
    int nNodes = nodeTable.size();
    V result = 0.0;
    V locResult = 0.0;

    for (int n = 0; n < nNodes; n++) {
        const auto &braNode = static_cast<const FunctionNode<D, T> &>(*nodeTable[n]);
        const MWNode<D, U> *mwNode = ket.findNode(braNode.getNodeIndex());
        if (mwNode == nullptr) continue;

        const auto &ketNode = static_cast<const FunctionNode<D, U> &>(*mwNode);
        if (braNode.isRootNode()) locResult += dot_scaling(braNode, ketNode);
        locResult += dot_wavelet(braNode, ketNode);
    }
    result += locResult;
    return result;
}

/**
 * @brief Absolute inner product proxy based on node norms.
 *
 * @tparam D  Spatial dimension.
 * @tparam T  Coefficient type.
 *
 * @param[in] bra   First input function.
 * @param[in] ket   Second input function.
 * @param[in] exact If true, requires ket's grid to include bra's grid and
 *                  uses absolute coefficients per node. If false, uses an
 *                  approximate product of node norms and root-node norms.
 *
 * @returns Value proportional to the absolute inner product.
 *
 * @details
 * With exact = true, the routine converts to interpolating coefficients,
 * takes absolute values, and accumulates exact contributions node by node.
 * With exact = false, it avoids per-coefficient access and approximates the
 * product via node norms; disjoint functions yield zero.
 */
template <int D, typename T>
double node_norm_dot(FunctionTree<D, T> &bra, FunctionTree<D, T> &ket, bool exact) {
    if (bra.getMRA() != ket.getMRA()) MSG_ABORT("Incompatible MRA");

    double result = 0.0;
    int ncoef = bra.getKp1_d() * bra.getTDim();
    T valA[ncoef];
    T valB[ncoef];
    int nNodes = bra.getNEndNodes();

    for (int n = 0; n < nNodes; n++) {
        FunctionNode<D, T> &node = bra.getEndFuncNode(n);
        const NodeIndex<D> idx = node.getNodeIndex();
        if (exact) {
            FunctionNode<D, T> *mwNode = static_cast<FunctionNode<D, T> *>(ket.findNode(idx));
            if (mwNode == nullptr) MSG_ABORT("Trees must have same grid");
            node.getAbsCoefs(valA);
            mwNode->getAbsCoefs(valB);
            for (int i = 0; i < ncoef; i++) result += std::norm(valA[i] * valB[i]);
        } else {
            int rIdx = ket.getRootBox().getBoxIndex(idx);
            assert(rIdx >= 0);
            const MWNode<D, T> &root = ket.getRootBox().getNode(rIdx);
            result += std::sqrt(node.getSquareNorm()) * root.getNodeNorm(idx);
        }
    }

    return result;
}

// ---- Explicit instantiations ------------------------------------------------

template void
multiply<1, double>(double prec, FunctionTree<1, double> &out, double c, FunctionTree<1, double> &tree_a, FunctionTree<1, double> &tree_b, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void
multiply<2, double>(double prec, FunctionTree<2, double> &out, double c, FunctionTree<2, double> &tree_a, FunctionTree<2, double> &tree_b, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void
multiply<3, double>(double prec, FunctionTree<3, double> &out, double c, FunctionTree<3, double> &tree_a, FunctionTree<3, double> &tree_b, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply<1, double>(double prec, FunctionTree<1, double> &out, FunctionTreeVector<1, double> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply<2, double>(double prec, FunctionTree<2, double> &out, FunctionTreeVector<2, double> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply<3, double>(double prec, FunctionTree<3, double> &out, FunctionTreeVector<3, double> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply<1, double>(double prec, FunctionTree<1, double> &out, std::vector<FunctionTree<1, double> *> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply<2, double>(double prec, FunctionTree<2, double> &out, std::vector<FunctionTree<2, double> *> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply<3, double>(double prec, FunctionTree<3, double> &out, std::vector<FunctionTree<3, double> *> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void power<1, double>(double prec, FunctionTree<1, double> &out, FunctionTree<1, double> &tree, double pow, int maxIter, bool absPrec);
template void power<2, double>(double prec, FunctionTree<2, double> &out, FunctionTree<2, double> &tree, double pow, int maxIter, bool absPrec);
template void power<3, double>(double prec, FunctionTree<3, double> &out, FunctionTree<3, double> &tree, double pow, int maxIter, bool absPrec);
template void square<1, double>(double prec, FunctionTree<1, double> &out, FunctionTree<1, double> &tree, int maxIter, bool absPrec, bool conjugate);
template void square<2, double>(double prec, FunctionTree<2, double> &out, FunctionTree<2, double> &tree, int maxIter, bool absPrec, bool conjugate);
template void square<3, double>(double prec, FunctionTree<3, double> &out, FunctionTree<3, double> &tree, int maxIter, bool absPrec, bool conjugate);
template void dot<1, double>(double prec, FunctionTree<1, double> &out, FunctionTreeVector<1, double> &inp_a, FunctionTreeVector<1, double> &inp_b, int maxIter, bool absPrec);
template void dot<2, double>(double prec, FunctionTree<2, double> &out, FunctionTreeVector<2, double> &inp_a, FunctionTreeVector<2, double> &inp_b, int maxIter, bool absPrec);
template void dot<3, double>(double prec, FunctionTree<3, double> &out, FunctionTreeVector<3, double> &inp_a, FunctionTreeVector<3, double> &inp_b, int maxIter, bool absPrec);
template double node_norm_dot<1, double>(FunctionTree<1, double> &bra, FunctionTree<1, double> &ket, bool exact);
template double node_norm_dot<2, double>(FunctionTree<2, double> &bra, FunctionTree<2, double> &ket, bool exact);
template double node_norm_dot<3, double>(FunctionTree<3, double> &bra, FunctionTree<3, double> &ket, bool exact);

template void multiply<1, ComplexDouble>(double prec,
                                         FunctionTree<1, ComplexDouble> &out,
                                         ComplexDouble c,
                                         FunctionTree<1, ComplexDouble> &tree_a,
                                         FunctionTree<1, ComplexDouble> &tree_b,
                                         int maxIter,
                                         bool absPrec,
                                         bool useMaxNorms,
                                         bool conjugate);
template void multiply<2, ComplexDouble>(double prec,
                                         FunctionTree<2, ComplexDouble> &out,
                                         ComplexDouble c,
                                         FunctionTree<2, ComplexDouble> &tree_a,
                                         FunctionTree<2, ComplexDouble> &tree_b,
                                         int maxIter,
                                         bool absPrec,
                                         bool useMaxNorms,
                                         bool conjugate);
template void multiply<3, ComplexDouble>(double prec,
                                         FunctionTree<3, ComplexDouble> &out,
                                         ComplexDouble c,
                                         FunctionTree<3, ComplexDouble> &tree_a,
                                         FunctionTree<3, ComplexDouble> &tree_b,
                                         int maxIter,
                                         bool absPrec,
                                         bool useMaxNorms,
                                         bool conjugate);
template void multiply<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, FunctionTreeVector<1, ComplexDouble> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, FunctionTreeVector<2, ComplexDouble> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, FunctionTreeVector<3, ComplexDouble> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void
multiply<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, std::vector<FunctionTree<1, ComplexDouble> *> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void
multiply<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, std::vector<FunctionTree<2, ComplexDouble> *> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void
multiply<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, std::vector<FunctionTree<3, ComplexDouble> *> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate);
template void power<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &tree, double pow, int maxIter, bool absPrec);
template void power<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &tree, double pow, int maxIter, bool absPrec);
template void power<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &tree, double pow, int maxIter, bool absPrec);
template void square<1, ComplexDouble>(double prec, FunctionTree<1, ComplexDouble> &out, FunctionTree<1, ComplexDouble> &tree, int maxIter, bool absPrec, bool conjugate);
template void square<2, ComplexDouble>(double prec, FunctionTree<2, ComplexDouble> &out, FunctionTree<2, ComplexDouble> &tree, int maxIter, bool absPrec, bool conjugate);
template void square<3, ComplexDouble>(double prec, FunctionTree<3, ComplexDouble> &out, FunctionTree<3, ComplexDouble> &tree, int maxIter, bool absPrec, bool conjugate);
template void dot<1, ComplexDouble>(double prec,
                                    FunctionTree<1, ComplexDouble> &out,
                                    FunctionTreeVector<1, ComplexDouble> &inp_a,
                                    FunctionTreeVector<1, ComplexDouble> &inp_b,
                                    int maxIter,
                                    bool absPrec);
template void dot<2, ComplexDouble>(double prec,
                                    FunctionTree<2, ComplexDouble> &out,
                                    FunctionTreeVector<2, ComplexDouble> &inp_a,
                                    FunctionTreeVector<2, ComplexDouble> &inp_b,
                                    int maxIter,
                                    bool absPrec);
template void dot<3, ComplexDouble>(double prec,
                                    FunctionTree<3, ComplexDouble> &out,
                                    FunctionTreeVector<3, ComplexDouble> &inp_a,
                                    FunctionTreeVector<3, ComplexDouble> &inp_b,
                                    int maxIter,
                                    bool absPrec);

template double dot<1, double, double>(FunctionTree<1, double> &bra, FunctionTree<1, double> &ket);
template double dot<2, double, double>(FunctionTree<2, double> &bra, FunctionTree<2, double> &ket);
template double dot<3, double, double>(FunctionTree<3, double> &bra, FunctionTree<3, double> &ket);
template ComplexDouble dot<1, ComplexDouble, double>(FunctionTree<1, ComplexDouble> &bra, FunctionTree<1, double> &ket);
template ComplexDouble dot<2, ComplexDouble, double>(FunctionTree<2, ComplexDouble> &bra, FunctionTree<2, double> &ket);
template ComplexDouble dot<3, ComplexDouble, double>(FunctionTree<3, ComplexDouble> &bra, FunctionTree<3, double> &ket);
template ComplexDouble dot<1, double, ComplexDouble>(FunctionTree<1, double> &bra, FunctionTree<1, ComplexDouble> &ket);
template ComplexDouble dot<2, double, ComplexDouble>(FunctionTree<2, double> &bra, FunctionTree<2, ComplexDouble> &ket);
template ComplexDouble dot<3, double, ComplexDouble>(FunctionTree<3, double> &bra, FunctionTree<3, ComplexDouble> &ket);
template ComplexDouble dot<1, ComplexDouble, ComplexDouble>(FunctionTree<1, ComplexDouble> &bra, FunctionTree<1, ComplexDouble> &ket);
template ComplexDouble dot<2, ComplexDouble, ComplexDouble>(FunctionTree<2, ComplexDouble> &bra, FunctionTree<2, ComplexDouble> &ket);
template ComplexDouble dot<3, ComplexDouble, ComplexDouble>(FunctionTree<3, ComplexDouble> &bra, FunctionTree<3, ComplexDouble> &ket);

template double node_norm_dot<1, ComplexDouble>(FunctionTree<1, ComplexDouble> &bra, FunctionTree<1, ComplexDouble> &ket, bool exact);
template double node_norm_dot<2, ComplexDouble>(FunctionTree<2, ComplexDouble> &bra, FunctionTree<2, ComplexDouble> &ket, bool exact);
template double node_norm_dot<3, ComplexDouble>(FunctionTree<3, ComplexDouble> &bra, FunctionTree<3, ComplexDouble> &ket, bool exact);

} // namespace mrcpp
