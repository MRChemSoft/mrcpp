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

#include <algorithm>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrcpp {

/* ----------------------------- Time evolution ----------------------------- */

void TimeEvolution_CrossCorrelationCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();
    int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol:
            MSG_ERROR("Invalid scaling type");
            break;
        case Legendre:
            applyCcc(node);
            break;
        default:
            MSG_ERROR("Invalid scaling type");
            break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

void TimeEvolution_CrossCorrelationCalculator::applyCcc(MWNode<2> &node) {
    const int t_dim  = node.getTDim();   // 4
    const int kp1_d  = node.getKp1_d();  // (k+1)^2
    const int kmax   = node.getOrder();

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);

    // J at level (scale + 1)
    auto &Jlvl = *this->J_power_inetgarls[node.getScale() + 1];
    const size_t span = Jlvl.size(); // number of l-shifts available (2N-1)

    const NodeIndex<2> &idx = node.getNodeIndex();

    for (int i = 0; i < t_dim; ++i) {
        NodeIndex<2> l = idx.child(i);
        int l_b = l[1] - l[0];

        // wrap to [0, span)
        size_t l_idx = (l_b >= 0) ? static_cast<size_t>(l_b)
                                  : static_cast<size_t>(span + l_b);
        if (l_idx >= span) continue;

        const auto &Jvec = Jlvl[l_idx];          // std::vector<std::complex<double>>
        const int Jlen   = static_cast<int>(Jvec.size());
        const int K_avail = static_cast<int>(this->cross_correlation->Matrix.size()) - 1;

        int pos = 0;
        for (int p = 0; p <= kmax; ++p) {
            for (int j = 0; j <= kmax; ++j, ++pos) {
                // m = 2k + p + j  < Jlen  => k ≤ (Jlen-1 - p - j)/2
                int K_by_J = (Jlen - 1 - p - j) / 2;
                int Kuse   = std::min(K_by_J, K_avail);
                if (Kuse < 0) continue;

                for (int k = 0; k <= Kuse; ++k) {
                    const int m = 2 * k + p + j;
                    double J = this->imaginary ? Jvec[m].imag() : Jvec[m].real();
                    vec_o.segment(i * kp1_d, kp1_d)(pos) +=
                        J * this->cross_correlation->Matrix[k](p, j);
                }
            }
        }
    }

    double *coefs = node.getCoefs();
    for (int n = 0; n < t_dim * kp1_d; ++n) coefs[n] = vec_o(n);
}

/* ------------------------------ Derivative -------------------------------- */

void DerivativeCrossCorrelationCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();
    int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol:
            MSG_ERROR("Invalid scaling type");
            break;
        case Legendre:
            applyCcc(node);
            break;
        default:
            MSG_ERROR("Invalid scaling type");
            break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

void DerivativeCrossCorrelationCalculator::applyCcc(MWNode<2> &node) {
    const int t_dim  = node.getTDim();   // 4
    const int kp1_d  = node.getKp1_d();  // (k+1)^2
    const int kmax   = node.getOrder();

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);

    auto &Jlvl = *this->J_power_inetgarls[node.getScale() + 1]; // DerivativePowerIntegrals
    const size_t span = Jlvl.size();

    const NodeIndex<2> &idx = node.getNodeIndex();

    for (int i = 0; i < t_dim; ++i) {
        NodeIndex<2> l = idx.child(i);
        int l_b = l[1] - l[0];

        size_t l_idx = (l_b >= 0) ? static_cast<size_t>(l_b)
                                  : static_cast<size_t>(span + l_b);
        if (l_idx >= span) continue;

        const auto &Jvec = Jlvl[l_idx];     // std::vector<double>, indexed by m
        const int Jlen   = static_cast<int>(Jvec.size());
        const int K_avail = static_cast<int>(this->cross_correlation->Matrix.size()) - 1;

        int pos = 0;
        for (int p = 0; p <= kmax; ++p) {
            for (int j = 0; j <= kmax; ++j, ++pos) {
                // m = 2k + 1 + p + j < Jlen => k ≤ (Jlen - 2 - p - j) / 2
                int K_by_J = (Jlen - 2 - p - j) / 2;
                int Kuse   = std::min(K_by_J, K_avail);
                if (Kuse < 0) continue;

                for (int k = 0; k <= Kuse; ++k) {
                    const int m = 2 * k + 1 + p + j; // odd index for derivative
                    double J = Jvec[m];
                    vec_o.segment(i * kp1_d, kp1_d)(pos) +=
                        J * this->cross_correlation->Matrix[k](p, j);
                }
            }
        }
    }

    double *coefs = node.getCoefs();
    for (int n = 0; n < t_dim * kp1_d; ++n) coefs[n] = vec_o(n);
}

} // namespace mrcpp