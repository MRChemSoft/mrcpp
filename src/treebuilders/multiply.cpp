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

template <int D, typename T>
void multiply(double prec, FunctionTree<D, T> &out, T c, FunctionTree<D, T> &inp_a, FunctionTree<D, T> &inp_b, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate) {
    FunctionTreeVector<D, T> tmp_vec;
    tmp_vec.push_back({c, &inp_a});
    tmp_vec.push_back({1.0, &inp_b});
    multiply(prec, out, tmp_vec, maxIter, absPrec, useMaxNorms, conjugate);
}

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

template <int D, typename T>
void multiply(double prec, FunctionTree<D, T> &out, std::vector<FunctionTree<D, T> *> &inp, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate) {
    FunctionTreeVector<D, T> inp_vec;
    for (auto &t : inp) inp_vec.push_back({1.0, t});
    multiply(prec, out, inp_vec, maxIter, absPrec, useMaxNorms, conjugate);
}

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
