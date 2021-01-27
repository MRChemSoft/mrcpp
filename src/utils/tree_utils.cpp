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

#include "tree_utils.h"

#include <algorithm>
#include <cmath>

#include "MRCPP/constants.h"

#include "math_utils.h"
#include "Printer.h"

#include "trees/HilbertIterator.h"
#include "trees/MWTree.h"
#include "trees/MWNode.h"


namespace mrcpp {

/** Calculate the threshold for the wavelet norm.
*
* Calculates the threshold that has to be met in the wavelet norm in order to
* guarantee the precision in the function representation. Depends on the
* square norm of the function and the requested relative accuracy. */
template <int D> bool tree_utils::split_check(const MWNode<D> &node, double prec, double split_fac, bool abs_prec) {
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

/** Traverse tree along the Hilbert path and find nodes of any rankId.
 * Returns one nodeVector for the whole tree. GenNodes disregarded. */
template <int D> void tree_utils::make_node_table(MWTree<D> &tree, MWNodeVector<D> &table) {
    HilbertIterator<D> it(&tree);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        table.push_back(&node);
    }
}

/** Traverse tree along the Hilbert path and find nodes of any rankId.
 * Returns one nodeVector per scale. GenNodes disregarded. */
template <int D> void tree_utils::make_node_table(MWTree<D> &tree, std::vector<MWNodeVector<D>> &table) {
    HilbertIterator<D> it(&tree);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        int depth = node.getDepth();
        // Add one more element
        if (depth + 1 > table.size()) table.push_back(MWNodeVector<D>());
        table[depth].push_back(&node);
    }
}

/** Make children scaling coefficients from parent
 * Other node info are not used/set
 * coeff_in are not modified.
 * The output is written directly into the 8 children scaling coefficients.
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
template <int D> void tree_utils::mw_transform(const MWTree<D> &tree,
                                               double *coeff_in,
                                               double *coeff_out,
                                               bool readOnlyScaling,
                                               int stride,
                                               bool b_overwrite) {
    int operation = Reconstruction;
    int kp1 = tree.getKp1();
    int kp1_d = tree.getKp1_d();
    int tDim = (1 << D);
    int kp1_dm1 = math_utils::ipow(kp1, D - 1);
    const MWFilter &filter = tree.getMRA().getFilter();
    double overwrite = 0.0;
    double tmpcoeff[kp1_d * tDim];
    double tmpcoeff2[kp1_d * tDim];
    int ftlim = tDim;
    int ftlim2 = tDim;
    int ftlim3 = tDim;
    if (readOnlyScaling) {
        ftlim = 1;
        ftlim2 = 2;
        ftlim3 = 4;
        // NB: Careful: tmpcoeff tmpcoeff2 are not initialized to zero
        // must not read these unitialized values!
    }

    overwrite = 0.0;
    int i = 0;
    int mask = 1;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = tmpcoeff + gt * kp1_d;
        for (int ft = 0; ft < ftlim; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
                double *in = coeff_in + ft * kp1_d;
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
        mask = 2; // 1 << i;
        for (int gt = 0; gt < tDim; gt++) {
            double *out = tmpcoeff2 + gt * kp1_d;
            for (int ft = 0; ft < ftlim2; ft++) {
                // Operate in direction i only if the bits along other
                // directions are identical. The bit of the direction we
                // operate on determines the appropriate filter/operator
                if ((gt | mask) == (ft | mask)) {
                    double *in = tmpcoeff + ft * kp1_d;
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
        mask = 4; // 1 << i;
        for (int gt = 0; gt < tDim; gt++) {
            double *out = coeff_out + gt * stride; // write right into children
            for (int ft = 0; ft < ftlim3; ft++) {
                // Operate in direction i only if the bits along other
                // directions are identical. The bit of the direction we
                // operate on determines the appropriate filter/operator
                if ((gt | mask) == (ft | mask)) {
                    double *in = tmpcoeff2 + ft * kp1_d;
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
        double *out;
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
template <int D> void tree_utils::mw_transform_back(MWTree<D> &tree, double *coeff_in, double *coeff_out, int stride) {
    NOT_IMPLEMENTED_ABORT;
}

/** Make parent from children scaling coefficients
 * Other node info are not used/set
 * coeff_in are not modified.
 * The output is read directly from the 8 children scaling coefficients.
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
template <> void tree_utils::mw_transform_back<3>(MWTree<3> &tree, double *coeff_in, double *coeff_out, int stride) {
    int operation = Compression;
    int kp1 = tree.getKp1();
    int kp1_d = tree.getKp1_d();
    int tDim = 8;
    int kp1_dm1 = math_utils::ipow(kp1, 2);
    const MWFilter &filter = tree.getMRA().getFilter();
    double overwrite = 0.0;
    double tmpcoeff[kp1_d * tDim];

    int ftlim = tDim;
    int ftlim2 = tDim;
    int ftlim3 = tDim;

    int i = 0;
    int mask = 1;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = coeff_out + gt * kp1_d;
        for (int ft = 0; ft < ftlim; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
                double *in = coeff_in + ft * stride;
                int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    i++;
    mask = 2; // 1 << i;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = tmpcoeff + gt * kp1_d;
        for (int ft = 0; ft < ftlim2; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
                double *in = coeff_out + ft * kp1_d;
                int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    i++;
    mask = 4; // 1 << i;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = coeff_out + gt * kp1_d;
        // double *out = coeff_out + gt * N_coeff;
        for (int ft = 0; ft < ftlim3; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
                double *in = tmpcoeff + ft * kp1_d;
                int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                const Eigen::MatrixXd &oper = filter.getSubFilter(filter_index, operation);

                math_utils::apply_filter(out, in, oper, kp1, kp1_dm1, overwrite);
                overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
}

template bool tree_utils::split_check<1>(const MWNode<1> &node, double prec, double split_fac, bool abs_prec);
template bool tree_utils::split_check<2>(const MWNode<2> &node, double prec, double split_fac, bool abs_prec);
template bool tree_utils::split_check<3>(const MWNode<3> &node, double prec, double split_fac, bool abs_prec);

template void tree_utils::make_node_table<1>(MWTree<1> &tree, MWNodeVector<1> &table);
template void tree_utils::make_node_table<2>(MWTree<2> &tree, MWNodeVector<2> &table);
template void tree_utils::make_node_table<3>(MWTree<3> &tree, MWNodeVector<3> &table);

template void tree_utils::make_node_table<1>(MWTree<1> &tree, std::vector<MWNodeVector<1>> &table);
template void tree_utils::make_node_table<2>(MWTree<2> &tree, std::vector<MWNodeVector<2>> &table);
template void tree_utils::make_node_table<3>(MWTree<3> &tree, std::vector<MWNodeVector<3>> &table);

template void tree_utils::mw_transform<1>(const MWTree<1> &tree, double *coeff_in, double *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);
template void tree_utils::mw_transform<2>(const MWTree<2> &tree, double *coeff_in, double *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);
template void tree_utils::mw_transform<3>(const MWTree<3> &tree, double *coeff_in, double *coeff_out, bool readOnlyScaling, int stride, bool b_overwrite);

template void tree_utils::mw_transform_back<1>(MWTree<1> &tree, double *coeff_in, double *coeff_out, int stride);
template void tree_utils::mw_transform_back<2>(MWTree<2> &tree, double *coeff_in, double *coeff_out, int stride);
template void tree_utils::mw_transform_back<3>(MWTree<3> &tree, double *coeff_in, double *coeff_out, int stride);

} // namespace mrcpp
