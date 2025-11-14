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

/*=====================  TimeEvolution_CrossCorrelationCalculator  =====================*/

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
    // Dimensions
    const int t_dim = node.getTDim();        // 4
    const int kp1_d = node.getKp1_d();       // (k+1)^2
    const int order = node.getOrder();

    // Safety checks to avoid segfaults
    if (!this->cross_correlation) {
        MSG_ERROR("cross_correlation pointer is null");
        return;
    }
    const int C_blocks = static_cast<int>(this->cross_correlation->Matrix.size());
    if (C_blocks <= 0) {
        MSG_ERROR("cross_correlation->Matrix is empty");
        return;
    }

    // Fetch the integrals for this scale (+1 as per your design)
    const int scale_key = node.getScale() + 1;
    auto itS = this->J_power_inetgarls.find(scale_key);
    if (itS == this->J_power_inetgarls.end() || !(itS->second)) {
        MSG_ERROR("Missing JpowerIntegrals for scale key " << scale_key);
        return;
    }
    auto &J_power_integrals_at_scale = *(itS->second);

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();

    for (int i = 0; i < t_dim; i++) {
        NodeIndex<2> l = idx.child(i);
        int l_b = l[1] - l[0];  // index into the per-shift table

        // Get the vector of complex J for this shift; operator[] handles negatives
        auto &J_vec = J_power_integrals_at_scale[l_b];
        const int Jsz = static_cast<int>(J_vec.size());

        int vec_o_segment_index = 0;
        for (int p = 0; p <= order; p++) {
            for (int j = 0; j <= order; j++) {
                // Bound on k from the J table: 2*k + p + j < Jsz
                int maxK_from_J = (Jsz - 1 - (p + j)) / 2; // floor
                if (maxK_from_J < 0) {
                    ++vec_o_segment_index;
                    continue; // nothing to add for this (p,j)
                }

                // Also bound by number of correlation blocks available
                int maxK = std::min(maxK_from_J, C_blocks - 1);

                // Accumulate
                for (int k = 0; k <= maxK; ++k) {
                    const int Jindex = 2 * k + p + j;
                    double Jval = this->imaginary ? J_vec[Jindex].imag() : J_vec[Jindex].real();

                    // Guard the matrix block shape (should be KxK); Eigen will assert if out-of-bounds
                    vec_o.segment(i * kp1_d, kp1_d)(vec_o_segment_index) +=
                        Jval * cross_correlation->Matrix[k](p, j);
                }
                ++vec_o_segment_index;
            }
        }
    }

    // Write coefficients back
    double *coefs = node.getCoefs();
    for (int i = 0; i < t_dim * kp1_d; i++) {
        coefs[i] = vec_o(i);
    }
}

/*=====================  DerivativeCrossCorrelationCalculator  =====================*/

void DerivativeCrossCorrelationCalculator::calcNode(MWNode<2> &node) {
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

void DerivativeCrossCorrelationCalculator::applyCcc(MWNode<2> &node) {
#pragma omp critical
{
    // Dimensions
    const int t_dim = node.getTDim();        // 4
    const int kp1_d = node.getKp1_d();       // (k+1)^2
    const int order = node.getOrder();

    // Safety checks to avoid segfaults
    if (!this->cross_correlation) {
        MSG_ERROR("cross_correlation pointer is null");
        return;
    }
    const int C_blocks = static_cast<int>(this->cross_correlation->Matrix.size());
    if (C_blocks <= 0) {
        MSG_ERROR("cross_correlation->Matrix is empty");
        return;
    }

    // Fetch the integrals for this scale (+1 as per your design)
    const int scale_key = node.getScale() + 1;
    auto itS = this->J_power_inetgarls.find(scale_key);
    if (itS == this->J_power_inetgarls.end() || !(itS->second)) {
        MSG_ERROR("Missing DerivativePowerIntegrals for scale key " << scale_key);
        return;
    }
    auto &J_power_integrals_at_scale = *(itS->second);

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();

    for (int i = 0; i < t_dim; i++) {
        NodeIndex<2> l = idx.child(i);
        int l_b = l[1] - l[0];  // shift index

        // For derivative case we have a vector<double> per shift
        auto &J_vec = J_power_integrals_at_scale[l_b];
        const int Jsz = static_cast<int>(J_vec.size());

        int vec_o_segment_index = 0;
        for (int p = 0; p <= order; p++) {
            for (int j = 0; j <= order; j++) {
                // Bound from the J table: 2*k + 1 + p + j < Jsz
                int maxK_from_J = (Jsz - 2 - (p + j)) / 2; // floor
                if (maxK_from_J < 0) {
                    ++vec_o_segment_index;
                    continue; // nothing to add for this (p,j)
                }

                // Also bound by number of correlation blocks available
                int maxK = std::min(maxK_from_J, C_blocks - 1);

                for (int k = 0; k <= maxK; ++k) {
                    const int Jindex = 2 * k + 1 + p + j;
                    const double Jval = J_vec[Jindex];

                    vec_o.segment(i * kp1_d, kp1_d)(vec_o_segment_index) +=
                        Jval * cross_correlation->Matrix[k](p, j);
                }
                ++vec_o_segment_index;
            }
        }
    }

    // Write coefficients back
    double *coefs = node.getCoefs();
    for (int i = 0; i < t_dim * kp1_d; i++) {
        coefs[i] = vec_o(i);
    }
} // omp critical
}

} // namespace mrcpp