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

#include "tree_utils.h"

#include <algorithm>
#include <cmath>

#include "MRCPP/constants.h"

#include "Printer.h"
#include "math_utils.h"

#include "trees/MWNode.h"
#include "trees/MWTree.h"
#include "trees/TreeIterator.h"

namespace mrcpp {

template <int D, typename T> bool tree_utils::split_check(const MWNode<D, T> &node, double prec, double split_fac, bool abs_prec) {
    bool split = false;
    if (prec > 0.0) {
        double t_norm = 1.0;
        double sq_norm = node.getMWTree().getSquareNorm();
        if (sq_norm > 0.0 and not abs_prec) t_norm = std::sqrt(sq_norm);

        double scale_fac = 1.0;
        if (split_fac > MachineZero) {
            double expo = 0.5 * split_fac * (node.getScale() + 1);
            scale_fac = std::pow(2.0, -expo);
        }

        double w_thrs = std::max(2.0 * MachinePrec, prec * t_norm * scale_fac);
        double w_norm = std::sqrt(node.getWaveletNorm());
        if (w_norm > w_thrs) split = true;
    }
    return split;
}

template <int D, typename T> void tree_utils::make_node_table(MWTree<D, T> &tree, MWNodeVector<D, T> &table) {
    TreeIterator<D, T> it(tree, TopDown, Hilbert);
    it.setReturnGenNodes(false);
    while (it.nextParent()) {
        MWNode<D, T> &node = it.getNode();
        if (node.getDepth() == 0) continue;
        table.push_back(&node);
    }
    it.init(tree);
    while (it.next()) {
        MWNode<D, T> &node = it.getNode();
        table.push_back(&node);
    }
}

template <int D, typename T> void tree_utils::make_node_table(MWTree<D, T> &tree, std::vector<MWNodeVector<D, T>> &table) {
    TreeIterator<D, T> it(tree, TopDown, Hilbert);
    it.setReturnGenNodes(false);
    while (it.nextParent()) {
        MWNode<D, T> &node = it.getNode();
        if (node.getDepth() == 0) continue;
        int depth = node.getDepth() + tree.getNNegScales();
        if (depth + 1 > table.size()) table.push_back(MWNodeVector<D, T>());
        table[depth].push_back(&node);
    }
    it.init(tree);
    while (it.next()) {
        MWNode<D, T> &node = it.getNode();
        int depth = node.getDepth() + tree.getNNegScales();
        if (depth + 1 > table.size()) table.push_back(MWNodeVector<D, T>());
        table[depth].push_back(&node);
    }
}

template <int D, typename T> void tree_utils::mw_transform(const MWTree<D, T> &tree, T *coeff_in, T *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite) {
    int operation = Reconstruction;
    int kp1 = tree.getKp1();
    int kp1_d = tree.getKp1_d();
    int tDim = (1 << D);
    int kp1_dm1 = math_utils::ipow(kp1, D - 1);
    const MWFilter &filter = tree.getMRA().getFilter();
    double overwrite = 0.0;
    T tmpcoeff[kp1_d * tDim];
    T tmpcoeff2[kp1_d * tDim];
    int ftlim = tDim;
    int ftlim2 = tDim;
    int ftlim3 = tDim;
    if (readOnlyScaling) {
        ftlim = 1;
        ftlim2 = 2;
        ftlim3 = 4;
    }

    overwrite = 0.0;
    int i = 0;
    int mask = 1;
    for (int gt = 0; gt < tDim; gt++) {
        T *out = tmpcoeff + gt * kp1_d;
        for (int ft = 0; ft < ftlim; ft++) {
            if ((gt | mask) == (ft | mask)) {
                T *in = coeff_in + ft * kp1_d;
                int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    if (D > 1) {
        i++;
        mask = 2;
        for (int gt = 0; gt < tDim; gt++) {
            T *out = tmpcoeff2 + gt * kp1_d;
            for (int ft = 0; ft < ftlim2; ft++) {
                if ((gt | mask) == (ft | mask)) {
                    T *in = tmpcoeff + ft * kp1_d;
                    int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                    const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                    math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                    overwrite = 1.0;
                }
            }
            overwrite = 0.0;
        }
    }
    if (D > 2) {
        overwrite = 1.0;
        if (b_overwrite) overwrite = 0.0;
        i++;
        mask = 4;
        for (int gt = 0; gt < tDim; gt++) {
            T *out = coeff_out + gt * stride;
            for (int ft = 0; ft < ftlim3; ft++) {
                if ((gt | mask) == (ft | mask)) {
                    T *in = tmpcoeff2 + ft * kp1_d;
                    int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                    const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                    math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                    overwrite = 1.0;
                }
            }
            overwrite = 1.0;
            if (b_overwrite) overwrite = 0.0;
        }
    }

    if (D > 3) MSG_ABORT("D>3 NOT IMPLEMENTED for S_mwtransform");

    if (D < 3) {
        T *out;
        if (D == 1) out = tmpcoeff;
        if (D == 2) out = tmpcoeff2;
        if (b_overwrite) {
            for (int j = 0; j < tDim; j++) {
                for (int i = 0; i < kp1_d; i++) { coeff_out[i + j * stride] = out[i + j * kp1_d]; }
            }
        } else {
            for (int j = 0; j < tDim; j++) {
                for (int i = 0; i < kp1_d; i++) { coeff_out[i + j * stride] += out[i + j * kp1_d]; }
            }
        }
    }
}

// Specialized for D=3 below.
// template <int D, typename T> void tree_utils::mw_transform_back(MWTree<D, T> &tree, double *coeff_in, double *coeff_out, int stride) {
//    NOT_IMPLEMENTED_ABORT;
//}

template <typename T> void tree_utils::mw_transform_back(MWTree<3, T> &tree, T *coeff_in, T *coeff_out, int stride) {
    int operation = Compression;
    int kp1 = tree.getKp1();
    int kp1_d = tree.getKp1_d();
    int tDim = 8;
    int kp1_dm1 = math_utils::ipow(kp1, 2);
    const MWFilter &filter = tree.getMRA().getFilter();
    double overwrite = 0.0;
    T tmpcoeff[kp1_d * tDim];

    int ftlim = tDim;
    int ftlim2 = tDim;
    int ftlim3 = tDim;

    int i = 0;
    int mask = 1;
    for (int gt = 0; gt < tDim; gt++) {
        T *out = coeff_out + gt * kp1_d;
        for (int ft = 0; ft < ftlim; ft++) {
            if ((gt | mask) == (ft | mask)) {
                T *in = coeff_in + ft * stride;
                int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    i++;
    mask = 2;
    for (int gt = 0; gt < tDim; gt++) {
        T *out = tmpcoeff + gt * kp1_d;
        for (int ft = 0; ft < ftlim2; ft++) {
            if ((gt | mask) == (ft | mask)) {
                T *in = coeff_out + ft * kp1_d;
                int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    i++;
    mask = 4;
    for (int gt = 0; gt < tDim; gt++) {
        T *out = coeff_out + gt * kp1_d;
        for (int ft = 0; ft < ftlim3; ft++) {
            if ((gt | mask) == (ft | mask)) {
                T *in = tmpcoeff + ft * kp1_d;
                int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
}

template void tree_utils::make_node_table<1, double>(MWTree<1, double> &tree, MWNodeVector<1, double> &table);
template void tree_utils::make_node_table<2, double>(MWTree<2, double> &tree, MWNodeVector<2, double> &table);
template void tree_utils::make_node_table<3, double>(MWTree<3, double> &tree, MWNodeVector<3, double> &table);

template void tree_utils::make_node_table<1, double>(MWTree<1, double> &tree, std::vector<MWNodeVector<1, double>> &table);
template void tree_utils::make_node_table<2, double>(MWTree<2, double> &tree, std::vector<MWNodeVector<2, double>> &table);
template void tree_utils::make_node_table<3, double>(MWTree<3, double> &tree, std::vector<MWNodeVector<3, double>> &table);

template bool tree_utils::split_check<1, double>(const MWNode<1, double> &node, double prec, double split_fac, bool abs_prec);
template bool tree_utils::split_check<2, double>(const MWNode<2, double> &node, double prec, double split_fac, bool abs_prec);
template bool tree_utils::split_check<3, double>(const MWNode<3, double> &node, double prec, double split_fac, bool abs_prec);

template void tree_utils::mw_transform<1, double>(const MWTree<1, double> &tree, double *coeff_in, double *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);
template void tree_utils::mw_transform<2, double>(const MWTree<2, double> &tree, double *coeff_in, double *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);
template void tree_utils::mw_transform<3, double>(const MWTree<3, double> &tree, double *coeff_in, double *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);

// template void tree_utils::mw_transform_back<1, double>(MWTree<1, double> &tree, double *coeff_in, double *coeff_out, int stride);
// template void tree_utils::mw_transform_back<2, double>(MWTree<2, double> &tree, double *coeff_in, double *coeff_out, int stride);
template void tree_utils::mw_transform_back<double>(MWTree<3, double> &tree, double *coeff_in, double *coeff_out, int stride);

template void tree_utils::make_node_table<1, ComplexDouble>(MWTree<1, ComplexDouble> &tree, MWNodeVector<1, ComplexDouble> &table);
template void tree_utils::make_node_table<2, ComplexDouble>(MWTree<2, ComplexDouble> &tree, MWNodeVector<2, ComplexDouble> &table);
template void tree_utils::make_node_table<3, ComplexDouble>(MWTree<3, ComplexDouble> &tree, MWNodeVector<3, ComplexDouble> &table);

template void tree_utils::make_node_table<1, ComplexDouble>(MWTree<1, ComplexDouble> &tree, std::vector<MWNodeVector<1, ComplexDouble>> &table);
template void tree_utils::make_node_table<2, ComplexDouble>(MWTree<2, ComplexDouble> &tree, std::vector<MWNodeVector<2, ComplexDouble>> &table);
template void tree_utils::make_node_table<3, ComplexDouble>(MWTree<3, ComplexDouble> &tree, std::vector<MWNodeVector<3, ComplexDouble>> &table);

template bool tree_utils::split_check<1, ComplexDouble>(const MWNode<1, ComplexDouble> &node, double prec, double split_fac, bool abs_prec);
template bool tree_utils::split_check<2, ComplexDouble>(const MWNode<2, ComplexDouble> &node, double prec, double split_fac, bool abs_prec);
template bool tree_utils::split_check<3, ComplexDouble>(const MWNode<3, ComplexDouble> &node, double prec, double split_fac, bool abs_prec);

template void tree_utils::mw_transform<1, ComplexDouble>(const MWTree<1, ComplexDouble> &tree, ComplexDouble *coeff_in, ComplexDouble *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);
template void tree_utils::mw_transform<2, ComplexDouble>(const MWTree<2, ComplexDouble> &tree, ComplexDouble *coeff_in, ComplexDouble *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);
template void tree_utils::mw_transform<3, ComplexDouble>(const MWTree<3, ComplexDouble> &tree, ComplexDouble *coeff_in, ComplexDouble *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);

// template void tree_utils::mw_transform_back<1, ComplexDouble>(MWTree<1, ComplexDouble &tree, ComplexDouble *coeff_in, ComplexDouble *coeff_out, int stride);
// template void tree_utils::mw_transform_back<2, ComplexDouble>(MWTree<2, ComplexDouble &tree, ComplexDouble *coeff_in, ComplexDouble *coeff_out, int stride);
template void tree_utils::mw_transform_back<ComplexDouble>(MWTree<3, ComplexDouble> &tree, ComplexDouble *coeff_in, ComplexDouble *coeff_out, int stride);

} // namespace mrcpp
