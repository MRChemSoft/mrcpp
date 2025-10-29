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
 * @file TimeEvolution_CrossCorrelationCalculator.cpp
 * @brief Compute node-local coefficients for the time-evolution cross–correlation
 *        contribution on a 2D multiwavelet node.
 *
 * @details
 * This calculator assembles, for each 2D node, the coefficients produced by a
 * cross–correlation kernel combined with a set of precomputed power integrals
 * \( J_m \). Conceptually, for each child-translation \( i \in \{0,\dots,t\_dim-1\} \)
 * and for each local polynomial pair \( (p,j) \) (with \(0\le p,j \le k\)),
 * we accumulate
 *
 * \f[
 *   \mathrm{vec\_o}_{i,(p,j)}
 *   \;+=\;
 *   \sum_{k=0}^{K(i,p,j)}
 *     J_{2k+p+j}^{(l_b)}
 *     \; \cdot \;
 *     \mathrm{CC}[k](p,j) ,
 * \f]
 *
 * where:
 *  - \(k\) is the polynomial order (node order),
 *  - \(l_b\) is the child offset along the 1D index difference (second minus first),
 *  - \(\mathrm{CC}[k](p,j)\) are entries of the cross–correlation matrices
 *    (one per \(k\)), and
 *  - \(J_{m}^{(l_b)}\) are power integrals looked up from
 *    `J_power_inetgarls[scale+1][l_b][m]` (note: member name “inetgarls” is kept as-is).
 *
 * The result vector `vec_o` of length `t_dim * kp1_d` (with `t_dim = 4` in 2D and
 * `kp1_d = (k+1)^2`) is written into the node coefficient buffer. The node is then
 * compressed (`mwTransform(Compression)`), marked as having coefficients, and its
 * norms are updated.
 *
 * @note
 *  - Only the Legendre scaling basis is currently supported here; Interpol is rejected.
 *  - The member flag `imaginary` selects whether the imaginary or real parts of the
 *    \(J\)-integrals are used.
 *  - The code assumes the cross–correlation matrices have been pre-populated in
 *    `cross_correlation->Matrix[k]`, consistent with the node order.
 *  - No world-box rescaling is applied in this routine (values are directly assigned).
 *
 * @warning
 *  - The routine relies on external consistency:
 *      * `J_power_inetgarls[scale+1]` must exist for the node’s scale.
 *      * `J_power_inetgarls[...][l_b]` must cover all accessed indices `2*k+p+j`.
 *      * `cross_correlation->Matrix[k]` must be dimension-compatible with `(p,j)`.
 *  - If these invariants are violated, out-of-bounds access may occur upstream;
 *    the caller is responsible for preparing inputs correctly.
 */

#include "TimeEvolution_CrossCorrelationCalculator.h"
#include "trees/FunctionTree.h"
#include "trees/MWNode.h"
#include "utils/Printer.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrcpp {

/**
 * @brief Assemble time-evolution cross–correlation coefficients on a 2D node,
 *        then compress to MW form.
 *
 * @param[in,out] node The target multiwavelet node (D=2).
 *
 * @details
 * 1. Zero current coefficients.
 * 2. Dispatch based on scaling basis type:
 *    - **Legendre**: compute through #applyCcc.
 *    - **Interpol**: rejected (not implemented for this calculator).
 * 3. Compress (`mwTransform(Compression)`), mark coefficients present, and
 *    update node norms.
 */
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

/**
 * @brief Core assembly routine for Legendre scaling basis.
 *
 * @param[in,out] node The target 2D node.
 *
 * @details
 * Let `t_dim = node.getTDim()` (in 2D this is 4) and `kp1_d = (k+1)^2` with
 * `k = node.getOrder()`. For each child index `i` we compute its child index
 * difference `l_b = l[1] - l[0]` and accumulate
 *
 * \f[
 *   \mathrm{vec\_o}[i,(p,j)]
 *   \;+=\;
 *   \sum_{k=0}^{K}
 *     J_{2k+p+j}^{(l_b)} \cdot \mathrm{CC}[k](p,j),
 * \f]
 *
 * writing the final `vec_o` into the node coefficient buffer without further
 * rescaling in this routine. If `imaginary == true`, the imaginary parts of
 * the \(J\)-integrals are used; otherwise the real parts are used.
 *
 * @pre
 *  - `this->J_power_inetgarls[node.getScale() + 1]` is allocated and populated.
 *  - `cross_correlation->Matrix[k]` exists for all accessed `k` and is
 *    indexable at `(p,j)`, with `0 <= p,j <= node.getOrder()`.
 */
void TimeEvolution_CrossCorrelationCalculator::applyCcc(MWNode<2> &node) {
    // Node configuration
    int t_dim = node.getTDim();   // e.g. 4 in 2D
    int kp1_d = node.getKp1_d();  // (k + 1)^2

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();

    // Access precomputed J-power integrals for the node scale (+1 by convention).
    auto &J_power_inetgarls = *this->J_power_inetgarls[node.getScale() + 1];

    // Loop over children and local basis pairs (p, j)
    for (int i = 0; i < t_dim; i++) {
        NodeIndex<2> l = idx.child(i);
        int l_b = l[1] - l[0];

        int vec_o_segment_index = 0;
        for (int p = 0; p <= node.getOrder(); p++)
            for (int j = 0; j <= node.getOrder(); j++) {
                // Accumulate up to the largest admissible 2k+p+j supported by J_power_inetgarls[l_b]
                for (int k = 0; 2 * k + p + j < J_power_inetgarls[l_b].size(); k++) {
                    double J;
                    if (this->imaginary)
                        J = J_power_inetgarls[l_b][2 * k + p + j].imag();
                    else
                        J = J_power_inetgarls[l_b][2 * k + p + j].real();

                    // Note: Eigen reads matrices row-major from file by default in this setup;
                    //       hence the comment about transposition in the original code.
                    vec_o.segment(i * kp1_d, kp1_d)(vec_o_segment_index) +=
                        J * cross_correlation->Matrix[k](p, j);
                }
                vec_o_segment_index++;
            }
    }

    // Write assembled values into the node coefficient buffer (no additional scaling here).
    double *coefs = node.getCoefs();
    for (int i = 0; i < t_dim * kp1_d; i++) {
        coefs[i] = vec_o(i);
    }
}

} // namespace mrcpp
