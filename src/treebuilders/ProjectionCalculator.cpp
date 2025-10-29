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
 * @file ProjectionCalculator.cpp
 * @brief Compute scaling/wavelet coefficients by projecting an analytic (or
 *        otherwise representable) function onto the MW basis on a given node.
 *
 * @details
 * The projection proceeds by evaluating the input function at a set of
 * expanded child quadrature/collocation points associated with the node,
 * then transforming these samples into scaling coefficients and finally
 * compressing into the wavelet representation:
 *
 * 1. `node.getExpandedChildPts(exp_pts)` returns a D×N matrix of evaluation
 *    points in *local* node coordinates, where N equals `node.getNCoefs()`.
 * 2. Each point is rescaled by the per-dimension world-box scaling factors
 *    (`scaling_factor[d]`) so that the user function is evaluated in
 *    physical coordinates.
 * 3. The raw samples are written into the node coefficient buffer and
 *    converted to scaling coefficients via `cvTransform(Backward)`.
 * 4. `mwTransform(Compression)` moves the representation to compressed MW
 *    form (wavelets across scales, scaling on roots).
 * 5. Bookkeeping: mark coefficients present and update (square-)norms.
 *
 * The calculator is stateless across nodes; it assumes that the caller
 * (TreeBuilder) handles traversal and refinement decisions (via an adaptor).
 *
 * @note
 *  - The assertion `exp_pts.cols() == node.getNCoefs()` guards consistency
 *    between the quadrature layout and the node’s coefficient count.
 *  - `scaling_factor` is typically extracted from the world box and allows
 *    non-unit, per-axis domain scaling.
 *  - This implementation works for both real and complex coefficient types.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (`double` or `ComplexDouble`).
 */

#include "ProjectionCalculator.h"
#include "trees/MWNode.h"
#include <cassert>

using Eigen::MatrixXd;

namespace mrcpp {

/**
 * @brief Project a single node by sampling the input function on the node's
 *        expanded child grid and transforming samples into MW coefficients.
 *
 * @param[in,out] node The MW node whose coefficients are to be computed.
 *
 * @pre `node.getExpandedChildPts(exp_pts)` provides exactly `getNCoefs()`
 *      columns (asserted).
 * @post
 *  - Node coefficients represent the function in compressed MW form.
 *  - `node.setHasCoefs()` is set and node norms are updated.
 *
 * @implementation
 *  - Samples are taken at expanded child points, rescaled by
 *    `scaling_factor[d]`, and passed to `func->evalf`.
 *  - `cvTransform(Backward)` maps collocation values → scaling coefficients.
 *  - `mwTransform(Compression)` converts to MW compressed representation.
 */
template <int D, typename T>
void ProjectionCalculator<D, T>::calcNode(MWNode<D, T> &node) {
    MatrixXd exp_pts;
    node.getExpandedChildPts(exp_pts);

    assert(exp_pts.cols() == node.getNCoefs());

    Coord<D> r;
    T *coefs = node.getCoefs();
    for (int i = 0; i < node.getNCoefs(); i++) {
        for (int d = 0; d < D; d++) { r[d] = scaling_factor[d] * exp_pts(d, i); }
        coefs[i] = this->func->evalf(r);
    }

    node.cvTransform(Backward);     // collocation values -> scaling coefficients
    node.mwTransform(Compression);  // scaling/wavelet compression on the node
    node.setHasCoefs();             // mark that the node now owns valid coefs
    node.calcNorms();               // update norms for refinement/threshholding
}

/* --------------------------------------------------------------------------
 * Legacy (interpolating) variant
 *
 * The block below shows an older, somewhat faster interpolating approach
 * that assumes an interpolating scaling basis. It is kept as reference;
 * it performs quadrature using cached roots/weights and writes block
 * coefficients directly before compressing. Enable with care as it
 * assumes specific basis properties (Interpol).
 *
 * template<int D>
 * void ProjectionCalculator<D, T>::calcNode(MWNode<D, T> &node) { ... }
 * -------------------------------------------------------------------------- */

/* Old interpolating version, somewhat faster
template<int D>
void ProjectionCalculator<D, T>::calcNode(MWNode<D, T> &node) {
    const ScalingBasis &sf = node.getMWTree().getMRA().getScalingBasis();
    if (sf.getScalingType() != Interpol) {
        NOT_IMPLEMENTED_ABORT;
    }
    int quadratureOrder = sf.getQuadratureOrder();
    getQuadratureCache(qc);
    const VectorXd &pts = qc.getRoots(quadratureOrder);
    const VectorXd &wgts = qc.getWeights(quadratureOrder);

    double tmp_coefs[node.getNCoefs()];

    int scale = node.getScale();
    int kp1_d = node.getKp1_d();

    double scaleFactor = 1.0 / std::pow(2.0, scale + 1.0);
    double sqrtScaleFactor = std::sqrt(scaleFactor);
    double point[D];

    static int tDim = 1 << D;
    for (int cIdx = 0; cIdx < tDim; cIdx++) {
        NodeIndex<D> nIdx(node.getNodeIndex(), cIdx);
        const int *l = nIdx.getTranslation();

        int indexCounter[D];
        for (int i = 0; i < D; i++) {
            indexCounter[i] = 0;
        }

        for (int i = 0; i < kp1_d; i++) {
            double coef = 1.0;
            for (int j = 0; j < D; j++) {
                point[j] = scaleFactor * (pts(indexCounter[j]) + l[j]);
                coef *= std::sqrt(wgts(indexCounter[j])) * sqrtScaleFactor;
            }

            tmp_coefs[i] = coef * this->func->evalf(point);

            indexCounter[0]++;
            for (int j = 0; j < D - 1; j++) {
                if (indexCounter[j] == quadratureOrder) {
                    indexCounter[j] = 0;
                    indexCounter[j + 1]++;
                }
            }
        }
        node.setCoefBlock(cIdx, kp1_d, tmp_coefs);
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}
*/

/// Explicit template instantiations
template class ProjectionCalculator<1, double>;
template class ProjectionCalculator<2, double>;
template class ProjectionCalculator<3, double>;

template class ProjectionCalculator<1, ComplexDouble>;
template class ProjectionCalculator<2, ComplexDouble>;
template class ProjectionCalculator<3, ComplexDouble>;

} // namespace mrcpp
