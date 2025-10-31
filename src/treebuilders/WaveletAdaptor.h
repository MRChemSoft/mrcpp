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

#pragma once

#include "TreeAdaptor.h"
#include "utils/Printer.h"
#include "utils/tree_utils.h"

namespace mrcpp {

/**
 * @class WaveletAdaptor
 * @brief Refinement policy based on wavelet-norm error indicators.
 *
 * @tparam D Spatial dimension of the multiwavelet tree.
 * @tparam T Coefficient value type (e.g., double, ComplexDouble).
 *
 * @details
 * This adaptor decides whether a node should be *split* (refined) by
 * comparing its wavelet contribution against a precision target.
 * Internally it relies on @ref mrcpp::tree_utils::split_check, which
 * examines a node's (accumulated) wavelet norm relative to:
 *
 *  - a global precision @ref prec (optionally absolute via @ref absPrec),
 *  - a user-provided, index-dependent scaling @ref precFunc (defaults to 1),
 *  - an extra scale-dependent attenuation factor @ref splitFac
 *    (used to bias refinement with depth).
 *
 * If the threshold is exceeded, @ref splitNode requests refinement.
 *
 * @code{.cpp}
 * // Typical usage:
 * WaveletAdaptor<3, double> adapt(1e-6, 20);               // prec, maxScale
 * adapt.setPrecFunction([](const NodeIndex<3>&){ return 2.0; }); // tighten locally
 * TreeBuilder<3, double> builder;
 * builder.split(tree, adapt, false);                        // passCoefs = false
 * @endcode
 */
template <int D, typename T>
class WaveletAdaptor : public TreeAdaptor<D, T> {
public:
    /**
     * @brief Construct a wavelet-based adaptor.
     *
     * @param pr Global target precision (relative unless @p ap is true).
     * @param ms Maximum refinement scale (forwarded to @ref TreeAdaptor).
     * @param ap If true, interpret @p pr as an **absolute** tolerance;
     *           otherwise use a **relative** tolerance w.r.t. function norm.
     * @param sf Split-factor controlling depth bias (≥ 0). When > 0,
     *           the threshold is scaled by \f$2^{-0.5\,sf\,(s+1)}\f$ at scale s,
     *           encouraging deeper refinement only when warranted.
     */
    WaveletAdaptor(double pr, int ms, bool ap = false, double sf = 1.0)
            : TreeAdaptor<D, T>(ms)
            , absPrec(ap)
            , prec(pr)
            , splitFac(sf) {}

    /// @brief Virtual destructor.
    ~WaveletAdaptor() override = default;

    /**
     * @brief Provide a spatially varying precision multiplier.
     *
     * @param prec_func Function returning a factor (default 1.0) for
     *                  a given node index. The effective threshold becomes
     *                  `prec * prec_func(idx)` (plus depth scaling via @ref splitFac).
     *
     * @note Use this to tighten or relax refinement in specific regions,
     *       e.g. around features of interest.
     */
    void setPrecFunction(const std::function<double(const NodeIndex<D> &idx)> &prec_func) {
        this->precFunc = prec_func;
    }

protected:
    /// @brief If true, treat @ref prec as an absolute tolerance; otherwise relative.
    bool absPrec;

    /// @brief Base precision target used by the wavelet thresholding rule.
    double prec;

    /**
     * @brief Scale-dependent attenuation of the threshold.
     *
     * @details A positive value reduces the threshold with depth, making
     * refinement stricter at finer scales. Set to 0.0 to disable.
     */
    double splitFac;

    /**
     * @brief Per-node precision multiplier (defaults to identity).
     *
     * @details The effective threshold is `prec * precFunc(idx)` before
     * applying the depth-dependent @ref splitFac scaling.
     */
    std::function<double(const NodeIndex<D> &idx)> precFunc =
        [](const NodeIndex<D> & /*idx*/) { return 1.0; };

    /**
     * @brief Decide whether a node should be split.
     *
     * @param node The candidate node.
     * @return `true` if the node's wavelet norm exceeds the computed threshold.
     *
     * @details
     * Computes a local tolerance as:
     * \f[
     *   \tau = \text{prec} \times \text{precFunc}(\text{idx})
     * \f]
     * (relative to the function norm unless @ref absPrec is set),
     * then applies an additional depth-dependent factor governed by
     * @ref splitFac, and finally compares against the node's wavelet norm.
     */
    bool splitNode(const MWNode<D, T> &node) const override {
        auto precFac = this->precFunc(node.getNodeIndex()); // returns 1.0 by default
        return tree_utils::split_check(node, this->prec * precFac, this->splitFac, this->absPrec);
    }
};

} // namespace mrcpp