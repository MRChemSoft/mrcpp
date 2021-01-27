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

#include "CrossCorrelationCalculator.h"
#include "trees/FunctionTree.h"
#include "trees/MWNode.h"
#include "utils/Printer.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrcpp {

void CrossCorrelationCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();
    int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol: {
            getCrossCorrelationCache(Interpol, ccc);
            applyCcc(node, ccc);
            break;
        }
        case Legendre: {
            getCrossCorrelationCache(Legendre, ccc);
            applyCcc(node, ccc);
            break;
        }
        default:
            MSG_ERROR("Invalid scaling type");
            break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

template <int T> void CrossCorrelationCalculator::applyCcc(MWNode<2> &node, CrossCorrelationCache<T> &ccc) {
    const MatrixXd &lMat = ccc.getLMatrix(node.getOrder());
    const MatrixXd &rMat = ccc.getRMatrix(node.getOrder());

    int scale = node.getScale() + 1;
    int t_dim = node.getTDim();
    int kp1_d = node.getKp1_d();

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();
    for (int i = 0; i < t_dim; i++) {
        NodeIndex<2> l = idx.child(i);
        int l_a = l[1] - l[0] - 1;
        int l_b = l[1] - l[0];

        NodeIndex<1> idx_a(scale, {l_a});
        NodeIndex<1> idx_b(scale, {l_b});

        const MWNode<1> &node_a = this->kernel->getNode(idx_a);
        const MWNode<1> &node_b = this->kernel->getNode(idx_b);

        VectorXd vec_a;
        VectorXd vec_b;
        node_a.getCoefs(vec_a);
        node_b.getCoefs(vec_b);

        const VectorXd &seg_a = vec_a.segment(0, node_a.getKp1_d());
        const VectorXd &seg_b = vec_b.segment(0, node_b.getKp1_d());
        vec_o.segment(i * kp1_d, kp1_d) = (lMat * seg_a + rMat * seg_b);
    }
    double *coefs = node.getCoefs();
    double two_n = std::pow(2.0, -scale / 2.0);
    for (int i = 0; i < t_dim * kp1_d; i++) {
        auto sf = node.getMWTree().getMRA().getWorldBox().getScalingFactor(0);
        // This is only implemented for unifrom scaling factors
        // hence the zero TODO: make it work for non-unifrom scaling
        coefs[i] = std::sqrt(sf) * two_n * vec_o(i);
    }
}

} // namespace mrcpp
