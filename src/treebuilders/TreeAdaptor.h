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

/**
 * @file
 * @brief Generic adapter that decides whether tree nodes should be refined (split).
 *
 * @details
 * `TreeAdaptor` provides a lightweight, policy-style interface used by the tree
 * execution engine to determine which nodes of a @ref MWTree should be split
 * (refined). Concrete adaptors implement the decision rule in the protected
 * pure-virtual @ref splitNode method.
 *
 * Typical workflow:
 * 1. Construct a concrete adaptor (e.g., one that inspects norms, wavelet
 *    content, error estimates, etc.).
 * 2. Optionally set a maximum refinement scale with @ref setMaxScale.
 * 3. Call @ref splitNodeVector to examine an input list of nodes and append any
 *    newly-created children to an output list for further processing.
 *
 * The adaptor enforces two built-in guards in @ref splitNodeVector:
 * - **Branch nodes** (internal/structural nodes without coefficients) are
 *   skipped.
 * - Nodes deeper than the allowed scale threshold are skipped:
 *   `node.getScale() + 2 > maxScale`.
 *
 * @note The `+2` slack prevents overshooting the configured refinement ceiling
 *       when subsequent passes may still need headroom (implementation detail).
 */

#include "MRCPP/mrcpp_declarations.h"
#include "trees/MWNode.h"

namespace mrcpp {

/**
 * @class TreeAdaptor
 * @brief Abstract base class for node-refinement policies.
 *
 * @tparam D Spatial dimension of the tree.
 * @tparam T Coefficient value type (e.g., `double`, `ComplexDouble`).
 *
 * @details
 * Concrete adaptors derive from this class and implement @ref splitNode to
 * express the criterion that decides whether a given leaf node should be
 * refined. The base class owns the *maximum scale* guard and the helper that
 * performs splitting and collects children.
 */
template <int D, typename T>
class TreeAdaptor {
public:
    /**
     * @brief Construct with an initial maximum refinement scale.
     * @param ms Maximum scale (depth) allowed for refinement.
     *
     * Nodes at scales for which `node.getScale() + 2 > ms` will **not** be
     * split by @ref splitNodeVector, regardless of the policy decision.
     */
    explicit TreeAdaptor(int ms)
            : maxScale(ms) {}

    /// Virtual destructor (polymorphic base).
    virtual ~TreeAdaptor() = default;

    /**
     * @brief Change the maximum refinement scale.
     * @param ms New ceiling for refinement depth.
     *
     * @see maxScale
     */
    void setMaxScale(int ms) { this->maxScale = ms; }

    /**
     * @brief Apply the refinement policy to a batch of nodes and collect children.
     *
     * @param[out] out
     *   Vector that will receive pointers to the **newly created children**
     *   (across all nodes that are decided to be split).
     * @param[in] inp
     *   Vector of candidate nodes to be tested for splitting.
     *
     * @details
     * For each node in @p inp the routine:
     *  - skips the node if it is a **branch node** (see `MWNode::isBranchNode`);
     *  - enforces the scale guard `node.getScale() + 2 > maxScale`;
     *  - calls the policy @ref splitNode; if `true`, it creates the children
     *    (`node.createChildren(true)`) and appends them to @p out.
     *
     * @note Ownership of nodes remains with the tree; this function only pushes
     *       pointers to existing/newly created nodes into @p out.
     */
    void splitNodeVector(MWNodeVector<D, T> &out, MWNodeVector<D, T> &inp) const {
        for (int n = 0; n < inp.size(); n++) {
            MWNode<D, T> &node = *inp[n];

            // Skip structural nodes (no coefficients)
            if (node.isBranchNode()) continue;

            // Enforce maximum scale guard with a +2 safety margin
            if (node.getScale() + 2 > this->maxScale) continue;

            // Delegate the decision to the concrete adaptor
            if (splitNode(node)) {
                node.createChildren(true);
                for (int i = 0; i < node.getNChildren(); i++) {
                    out.push_back(&node.getMWChild(i));
                }
            }
        }
    }

protected:
    /**
     * @brief Maximum allowed refinement scale (depth) for newly created nodes.
     *
     * Nodes for which `node.getScale() + 2 > maxScale` are not considered for
     * splitting within @ref splitNodeVector.
     */
    int maxScale;

    /**
     * @brief Decide whether a given leaf node should be refined.
     *
     * @param node Candidate node (guaranteed non-branch and within scale guard).
     * @return `true` if the node must be split, `false` otherwise.
     *
     * @details
     * Derived classes implement this method to express an application-specific
     * refinement criterion (e.g., wavelet-norm threshold, operator bandwidth,
     * error estimator, etc.). This method must be **pure** (no side-effects)
     * with respect to the tree topology; @ref splitNodeVector performs the
     * actual splitting.
     */
    virtual bool splitNode(const MWNode<D, T> &node) const = 0;
};

} // namespace mrcpp