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
 * @file CopyAdaptor.cpp
 * @brief Tree adaptor that **copies** (follows) an existing grid structure,
 *        optionally widened by a user-specified bandwidth.
 *
 * @details
 * `mrcpp::CopyAdaptor` is a `TreeAdaptor` used with `TreeBuilder` to produce an
 * output function tree whose refinement pattern mirrors one or more **reference
 * trees**. It decides whether a node should be split solely by inspecting the
 * presence of the corresponding **children** (and their integer-neighbor shifts)
 * in the reference trees.
 *
 * This adaptor is typically used to:
 *  - replicate an input grid for **fixed-grid operations** (e.g., local
 *    derivative applies where no adaptivity is desired), and
 *  - **widen** the grid along selected directions to accommodate operators
 *    whose stencils reach into neighboring nodes (e.g., first/second derivative
 *    stencils). The widening is controlled by a per-dimension integer
 *    bandwidth \f$ \text{bandWidth}[d] \ge 0 \f$.
 *
 * ### Split criterion
 * For a candidate node `node` and each of its children `c` (in tensor-product
 * sense), the adaptor checks, for each dimension \f$d\in\{0,\dots,D-1\}\f$,
 * every integer shift \f$ \delta \in [-\text{bandWidth}[d],\text{bandWidth}[d]] \f$:
 *
 * 1. Form the child index `bwIdx = node.child(c)` and add the shift on the
 *    current dimension: `bwIdx[d] += δ`.
 * 2. If **any** reference `FunctionTree` contains that child index, the adaptor
 *    returns **true** (requesting the split).
 *
 * If no such child is found in any reference, the adaptor returns **false**.
 *
 * ### Notes
 * - This adaptor is **purely topological**; it does not inspect coefficients.
 * - If `bw == nullptr`, all bandwidths default to `0` (exact copy of the
 *   reference grid).
 * - The reference set can be a single tree or a vector of trees; the union of
 *   their reachable children (with bandwidth widening) drives the output grid.
 *
 * ### Example
 * @code
 * int bw[3] = {1, 0, 0};                // widen one node on each side in x
 * CopyAdaptor<3,double> pre(out_inp, maxScale, bw);
 * TreeBuilder<3,double> builder;
 * DefaultCalculator<3,double> calc;     // no-op; we only want to build the grid
 * builder.build(out, calc, pre, -1);    // fixed grid construction
 * @endcode
 */

#include "CopyAdaptor.h"

#include <tuple>

namespace mrcpp {

/**
 * @brief Construct a copy adaptor that follows a single reference tree.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type (`double` or `ComplexDouble`).
 * @param t   Reference function tree to follow.
 * @param ms  Maximum scale allowed for splitting (forwarded to `TreeAdaptor`).
 * @param bw  Optional pointer to an array of length `D` with per-dimension
 *            integer bandwidths. If `nullptr`, all bandwidths are set to `0`.
 *
 * @details
 * The adaptor will request a split whenever a corresponding child (possibly
 * shifted by up to `bw[d]` in each dimension) exists in the reference tree.
 */
template <int D, typename T>
CopyAdaptor<D, T>::CopyAdaptor(FunctionTree<D, T> &t, int ms, int *bw)
        : TreeAdaptor<D, T>(ms) {
    setBandWidth(bw);
    tree_vec.push_back(std::make_tuple(1.0, &t));
}

/**
 * @brief Construct a copy adaptor that follows the **union** of several trees.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param t   Vector of `(coef, tree*)` pairs; only the tree pointers matter for
 *            the splitting logic, the coefficients are ignored.
 * @param ms  Maximum scale allowed for splitting.
 * @param bw  Optional per-dimension bandwidth array. If `nullptr`, zeros.
 *
 * @details
 * A split is requested if **any** tree in `t` contains the candidate child
 * (within the bandwidth neighborhood).
 */
template <int D, typename T>
CopyAdaptor<D, T>::CopyAdaptor(FunctionTreeVector<D, T> &t, int ms, int *bw)
        : TreeAdaptor<D, T>(ms)
        , tree_vec(t) {
    setBandWidth(bw);
}

/**
 * @brief Set the per-dimension bandwidths used to widen the copied grid.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param bw Pointer to an integer array of length `D`. If `nullptr`, all
 *           bandwidths are set to `0`.
 *
 * @note Negative entries are treated as `0` by the caller contract; this
 *       function simply copies the values. The split loop ranges over
 *       `[-bandWidth[d], +bandWidth[d]]`.
 */
template <int D, typename T> void CopyAdaptor<D, T>::setBandWidth(int *bw) {
    for (int d = 0; d < D; d++) {
        if (bw != nullptr) {
            this->bandWidth[d] = bw[d];
        } else {
            this->bandWidth[d] = 0;
        }
    }
}

/**
 * @brief Decide whether a node should be split to mirror (and widen) a reference grid.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param node Candidate node in the output tree.
 * @return `true` if any reference tree contains a corresponding child
 *         (considering bandwidth shifts), `false` otherwise.
 *
 * @details
 * For each tensor child `c` of `node` and each dimension `d`, the method scans
 * integer offsets `bw ∈ [-bandWidth[d], +bandWidth[d]]`. The candidate child
 * index is formed as:
 *
 * @code
 * NodeIndex<D> bwIdx = idx.child(c);
 * bwIdx[d] += bw;
 * @endcode
 *
 * If any reference tree contains `bwIdx`, a split is requested immediately.
 * The search stops on the first positive hit.
 *
 * @complexity
 * \f$ \mathcal{O}\big(T \cdot C \cdot \prod_{d=0}^{D-1} (2\,\text{bandWidth}[d]+1)\big) \f$,
 * where `T` is the number of reference trees and `C` is the number of tensor
 * children per node.
 */
template <int D, typename T> bool CopyAdaptor<D, T>::splitNode(const MWNode<D, T> &node) const {
    const NodeIndex<D> &idx = node.getNodeIndex();
    for (int c = 0; c < node.getTDim(); c++) {
        for (int d = 0; d < D; d++) {
            for (int bw = -this->bandWidth[d]; bw <= this->bandWidth[d]; bw++) {
                NodeIndex<D> bwIdx = idx.child(c);
                bwIdx[d] += bw;
                for (int i = 0; i < this->tree_vec.size(); i++) {
                    const FunctionTree<D, T> &func_i = get_func(tree_vec, i);
                    const MWNode<D, T> *node_i = func_i.findNode(bwIdx);
                    if (node_i != nullptr) return true;
                }
            }
        }
    }
    return false;
}

// Explicit instantiations
template class CopyAdaptor<1, double>;
template class CopyAdaptor<2, double>;
template class CopyAdaptor<3, double>;

template class CopyAdaptor<1, ComplexDouble>;
template class CopyAdaptor<2, ComplexDouble>;
template class CopyAdaptor<3, ComplexDouble>;

} // namespace mrcpp
