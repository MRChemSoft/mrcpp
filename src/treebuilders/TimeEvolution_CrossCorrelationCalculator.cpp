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
    const int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol: MSG_ERROR("Invalid scaling type"); break;
        case Legendre: applyCcc(node); break;
        default:       MSG_ERROR("Invalid scaling type"); break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

void TimeEvolution_CrossCorrelationCalculator::applyCcc(MWNode<2> &node) {
    const int t_dim = node.getTDim();      // 4
    const int kp1_d = node.getKp1_d();     // (k+1)^2
    const int Kmax  = static_cast<int>(cross_correlation->Matrix.size());

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();

    auto &J_power = *this->J_power_inetgarls[node.getScale() + 1];

    for (int i = 0; i < t_dim; ++i) {
        NodeIndex<2> l = idx.child(i);
        const int l_b = l[1] - l[0];

        int vec_o_segment_index = 0;
        for (int p = 0; p <= node.getOrder(); ++p) {
            for (int j = 0; j <= node.getOrder(); ++j) {
                // safe in both dimensions: J and Matrix
                const int Jsize = static_cast<int>(J_power[l_b].size());
                for (int k = 0; k < Kmax; ++k) {
                    const int Jidx = 2 * k + p + j;
                    if (Jidx >= Jsize) break;
                    double Jval = this->imaginary ? J_power[l_b][Jidx].imag()
                                                  : J_power[l_b][Jidx].real();
                    vec_o.segment(i * kp1_d, kp1_d)(vec_o_segment_index)
                        += Jval * cross_correlation->Matrix[k](p, j);
                }
                ++vec_o_segment_index;
            }
        }
    }

    double *coefs = node.getCoefs();
    for (int i = 0; i < t_dim * kp1_d; ++i) coefs[i] = vec_o(i);
}

void DerivativeCrossCorrelationCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();
    const int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol: MSG_ERROR("Invalid scaling type"); break;
        case Legendre: applyCcc(node); break;
        default:       MSG_ERROR("Invalid scaling type"); break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

void DerivativeCrossCorrelationCalculator::applyCcc(MWNode<2> &node) {
    const int t_dim = node.getTDim();
    const int kp1_d = node.getKp1_d();
    const int Kmax  = static_cast<int>(cross_correlation->Matrix.size());

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();

    auto &J_power = *this->J_power_inetgarls[node.getScale() + 1];

    for (int i = 0; i < t_dim; ++i) {
        NodeIndex<2> l = idx.child(i);
        const int l_b = l[1] - l[0];

        int vec_o_segment_index = 0;
        for (int p = 0; p <= node.getOrder(); ++p) {
            for (int j = 0; j <= node.getOrder(); ++j) {
                const int Jsize = static_cast<int>(J_power[l_b].size());
                for (int k = 0; k < Kmax; ++k) {
                    const int Jidx = 2 * k + 1 + p + j;
                    if (Jidx >= Jsize) break;
                    const double Jval = J_power[l_b][Jidx];
                    vec_o.segment(i * kp1_d, kp1_d)(vec_o_segment_index)
                        += Jval * cross_correlation->Matrix[k](p, j);
                }
                ++vec_o_segment_index;
            }
        }
    }

    double *coefs = node.getCoefs();
    for (int i = 0; i < t_dim * kp1_d; ++i) coefs[i] = vec_o(i);
}

} // namespace mrcpp