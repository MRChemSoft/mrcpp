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
 * @brief Minimal adaptor that unconditionally (or never) splits nodes.
 *
 * @details
 * This adaptor implements the @ref TreeAdaptor interface with a constant
 * split policy: depending on a boolean flag passed at construction,
 * every node presented by the traversal will either:
 *  - be **always split** (when `split == true`), or
 *  - be **never split** (when `split == false`).
 *
 * This is useful for:
 *  - unit tests and benchmarks (forcing a fixed refinement pattern),
 *  - creating a uniform grid up to a maximum scale,
 *  - disabling refinement while still running a calculator over an existing grid.
 */

#include "TreeAdaptor.h"

namespace mrcpp {

/**
 * @class SplitAdaptor
 * @brief Constant split/no-split adaptor for tree refinement.
 *
 * @tparam D Spatial dimension of the tree.
 * @tparam T Scalar coefficient type (defaults to `double`).
 *
 * @details
 * The adaptor inherits the depth/scale controls from @ref TreeAdaptor (e.g.,
 * the *maximum scale* passed to the base constructor). The split decision
 * itself is independent of node content and simply mirrors the `split` flag.
 *
 * ### Refinement semantics
 * - If `split == true`, any node that the base class allows to be refined
 *   (e.g., below the maximum scale) will be marked for splitting.
 * - If `split == false`, no node will be split by this adaptor (even if below
 *   the maximum scale).
 */
template <int D, typename T = double>
class SplitAdaptor final : public TreeAdaptor<D, T> {
public:
    /**
     * @brief Construct a constant split adaptor.
     *
     * @param[in] ms Maximum scale (or equivalent depth limit) forwarded to
     *               @ref TreeAdaptor. Nodes at or beyond this scale will not be
     *               refined by the base logic regardless of @p sp.
     * @param[in] sp Split policy: `true` to always split (subject to base
     *               constraints), `false` to never split.
     */
    SplitAdaptor(int ms, bool sp)
            : TreeAdaptor<D, T>(ms)
            , split(sp) {}

private:
    /// Constant split policy applied to every visited node.
    bool split;

    /**
     * @brief Decide whether to split a node.
     *
     * @param[in] node Node under consideration (unused).
     * @return `true` if this adaptor is configured to split; otherwise `false`.
     *
     * @note The base class may still veto refinement (e.g., beyond max scale).
     */
    bool splitNode(const MWNode<D, T> &node) const override { return this->split; }
};

} // namespace mrcpp
