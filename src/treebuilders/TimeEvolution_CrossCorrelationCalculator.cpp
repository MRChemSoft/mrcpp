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

#include "TimeEvolution_CrossCorrelationCalculator.h"
#include "trees/FunctionTree.h"
#include "trees/MWNode.h"
#include "utils/Printer.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrcpp {

void TimeEvolution_CrossCorrelationCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();
    int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol: {
            MSG_ERROR("Invalid scaling type");
            break;
        }
        case Legendre: {
            applyCcc(node);
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

void TimeEvolution_CrossCorrelationCalculator::applyCcc(MWNode<2> &node) {
    int t_dim = node.getTDim();
    int kp1_d = node.getKp1_d();

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();

    auto &J_power_inetgarls = *this->J_power_inetgarls[node.getScale() + 1];

    for (int i = 0; i < t_dim; i++) {
        NodeIndex<2> l = idx.child(i);
        int l_b = l[1] - l[0];

        int vec_o_segment_index = 0;
        for (int p = 0; p <= node.getOrder(); p++)
            for (int j = 0; j <= node.getOrder(); j++) {
                for (int k = 0; 2 * k + p + j < J_power_inetgarls[l_b].size(); k++) {
                    double J;
                    if (this->imaginary)
                        J = J_power_inetgarls[l_b][2 * k + p + j].imag();
                    else
                        J = J_power_inetgarls[l_b][2 * k + p + j].real();

                    vec_o.segment(i * kp1_d, kp1_d)(vec_o_segment_index) +=
                        J * cross_correlation->Matrix[k](p, j);
                }
                vec_o_segment_index++;
            }
    }

    double *coefs = node.getCoefs();
    for (int i = 0; i < t_dim * kp1_d; i++) {
        coefs[i] = vec_o(i);
    }
}

} // namespace mrcpp
