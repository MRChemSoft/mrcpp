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
 * @file CrossCorrelationCalculator.cpp
 * @brief Assembly of 2D cross–correlation operator blocks from a 1D kernel,
 *        for Legendre and Interpolating scaling bases.
 *
 * @details
 * This module implements the node-wise assembly of a separable 2D operator that
 * represents the *cross–correlation* between adjacent 1D kernel segments.
 * Given a 1D kernel stored as an `MWTree<1>` (accessed via
 * `CrossCorrelationCalculator::kernel`), the calculator:
 *
 *  - selects the appropriate **precomputed** cross–correlation matrices
 *    \f$L\f$ and \f$R\f$ from a `CrossCorrelationCache` depending on the
 *    scaling basis (Legendre or Interpolating) and the local polynomial order,
 *  - extracts the relevant 1D kernel coefficient blocks at indices shifted by
 *    the child offset of the current 2D node,
 *  - forms the 2D block by the linear combination
 *    \f[
 *      \mathbf{v}_o^{(i)} \;=\; L \,\mathbf{v}_a \;+\; R \,\mathbf{v}_b,
 *    \f]
 *    where \f$\mathbf{v}_a\f$ and \f$\mathbf{v}_b\f$ are the 1D kernel
 *    coefficient segments corresponding to the left/right neighboring child
 *    positions induced by the current 2D node child \f$i\f$,
 *  - applies a scale factor \f$ 2^{-\,(\text{scale}+1)/2} \f$ and a global
 *    normalization factor derived from the world-box scaling to obtain the
 *    final coefficient block for the 2D operator node.
 *
 * The assembled coefficients are then compressed (wavelet transform in
 * `Compression` mode), marked as present, and their norms are computed for
 * downstream thresholding and application.
 *
 * ### Indexing convention
 * For a node with index \f$\ell = (\ell_0,\ell_1)\f$ and a specific tensor
 * child \f$i\f$, the child index `l = idx.child(i)` induces two 1D offsets
 * \f[
 *     \ell_a = \ell_1 - \ell_0 - 1,
 *     \qquad
 *     \ell_b = \ell_1 - \ell_0,
 * \f]
 * which select adjacent 1D kernel nodes at the next finer scale. These are
 * mapped to 1D node indices \f$(\text{scale}+1,\ell_a)\f$ and
 * \f$(\text{scale}+1,\ell_b)\f$.
 *
 * @note At the moment, only **uniform scaling factors** are supported; the code
 *       reads the scaling factor for dimension 0 and assumes it is uniform.
 */

#include "CrossCorrelationCalculator.h"
#include "trees/FunctionTree.h"
#include "trees/MWNode.h"
#include "utils/Printer.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrcpp {

/**
 * @brief Build the cross–correlation block for a 2D operator node.
 *
 * @param node Output 2D operator node to be filled (overwrites coefficients).
 *
 * @details
 * - Zeros existing coefficients.
 * - Detects the scaling basis of the underlying MRA (`Interpol` or `Legendre`).
 * - Retrieves the corresponding `CrossCorrelationCache` and dispatches to
 *   `applyCcc()` which performs the actual linear algebra using cached
 *   matrices \f$L,R\f$ and the 1D kernel tree referenced by this calculator.
 * - Applies compression (`mwTransform(Compression)`), marks coefficients as
 *   present, and computes norms (`calcNorms()`).
 *
 * @throws Emits an error if the scaling type is unsupported.
 */
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

/**
 * @brief Assemble one 2D node from cached cross–correlation matrices and a 1D kernel.
 *
 * @tparam T Tag of the scaling basis cache (`Interpol` or `Legendre`).
 * @param node Output 2D operator node to be filled.
 * @param ccc  Cross–correlation cache for the selected scaling basis.
 *
 * @details
 * Let \f$k\f$ denote the 1D kernel function tree pointed to by
 * `this->kernel`. For each tensor child \f$i\f$ of the current 2D node:
 *  - compute the child index \f$l = \text{idx.child}(i)\f$ where
 *    \f$\text{idx}\f$ is the node index,
 *  - form the adjacent 1D indices \f$\ell_a = l_1-l_0-1\f$ and
 *    \f$\ell_b = l_1-l_0\f$ at scale \f$s = \text{node.getScale()}+1\f$,
 *  - fetch the 1D coefficient vectors \f$\mathbf{v}_a,\mathbf{v}_b\f$ from
 *    the kernel tree at \f$(s,\ell_a)\f$ and \f$(s,\ell_b)\f$,
 *  - compute the 2D segment
 *    \f[
 *       \mathbf{v}_o^{(i)} \;=\; L \,\mathbf{v}_a \;+\; R \,\mathbf{v}_b,
 *    \f]
 *    where \f$L,R\f$ are read from the cache for the node order
 *    (`ccc.getLMatrix(node.getOrder())`, `ccc.getRMatrix(node.getOrder())`),
 *  - store \f$\mathbf{v}_o^{(i)}\f$ in the appropriate slot of the 2D node
 *    coefficient buffer after applying the normalization
 *    \f$ \sqrt{\text{scaling\_factor}}\, 2^{-s/2} \f$.
 *
 * The method writes directly into `node.getCoefs()` and does not allocate
 * intermediate node structures beyond temporary vectors.
 *
 * @note Only uniform world-box scaling factors are supported at present.
 */
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

        VectorXd vec_a, vec_b;
        node_a.getCoefs(vec_a);
        node_b.getCoefs(vec_b);

        const VectorXd &seg_a = vec_a.segment(0, node_a.getKp1_d());
        const VectorXd &seg_b = vec_b.segment(0, node_b.getKp1_d());
        vec_o.segment(i * kp1_d, kp1_d) = (lMat * seg_a + rMat * seg_b);
    }
    double *coefs = node.getCoefs();
    double two_n = std::pow(2.0, -scale / 2.0);
    for (int i = 0; i < t_dim * kp1_d; i++) {
        auto scaling_factor = node.getMWTree().getMRA().getWorldBox().getScalingFactor(0);
        // Implemented for uniform scaling factors (dimension 0). For non-uniform
        // scaling a per-dimension normalization would be required.
        coefs[i] = std::sqrt(scaling_factor) * two_n * vec_o(i);
    }
}

} // namespace mrcpp
