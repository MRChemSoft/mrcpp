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
 * @file NodeBox.h
 * @brief Container that associates a regular grid of boxes with pointers to MW nodes.
 *
 * @details
 * A NodeBox is a thin wrapper around @ref BoundingBox that, in addition to the
 * geometric information (bounds, scale, periodicity), keeps a dense array of
 * pointers to @ref MWNode objects—one slot per box at the underlying scale.
 * It is used by @ref MWTree to store and access the set of **root nodes**
 * at the world scale, and by other components whenever a compact mapping
 * from box indices to nodes is required.
 */

#pragma once

#include "BoundingBox.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @class NodeBox
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Scalar type of the associated @ref MWNode (e.g., double, ComplexDouble).
 *
 * @brief Bounding box with node-pointer storage.
 *
 * @details
 * The class allocates and owns a contiguous array of pointers, one per box
 * defined by the base @ref BoundingBox. Pointers are not owned by NodeBox
 * (ownership stays with the corresponding @ref MWTree allocator); NodeBox only
 * stores and clears them. The counter @ref nOccupied tracks how many slots
 * are non-null.
 */
template <int D, typename T> class NodeBox final : public BoundingBox<D> {
public:
    /**
     * @brief Construct a NodeBox from a lower-corner index and number of boxes.
     * @param idx Lower-corner @ref NodeIndex at the world scale.
     * @param nb  Number of boxes per dimension (defaults to all ones).
     *
     * @details The geometric information is taken from @p idx and @p nb.
     * Internal pointer storage is allocated and initialized to `nullptr`.
     */
    NodeBox(const NodeIndex<D> &idx, const std::array<int, D> &nb = {});

    /**
     * @brief Copy-construct from another NodeBox.
     * @param box Source NodeBox.
     *
     * @details Copies the underlying @ref BoundingBox state and recreates
     * pointer storage; node pointers themselves are copied (shallow).
     */
    NodeBox(const NodeBox<D, T> &box);

    /**
     * @brief Construct from a plain @ref BoundingBox.
     * @param box Geometric box to take as base.
     *
     * @details Creates an equivalent NodeBox and allocates empty pointer storage.
     */
    NodeBox(const BoundingBox<D> &box);

    /// Non-assignable (pointer storage is managed per-instance).
    NodeBox<D, T> &operator=(const NodeBox<D, T> &box) = delete;

    /// Destructor; releases the internal pointer array (not the nodes).
    ~NodeBox() override;

    /**
     * @brief Store a node pointer in slot @p idx.
     * @param idx  Linear box index in `[0, size())`.
     * @param node Address of the node pointer to store (double pointer).
     *
     * @details The stored value is `*node`. If it was previously `nullptr`
     * and the new value is non-null, @ref nOccupied is incremented. If it was
     * non-null and is reset to `nullptr`, @ref nOccupied is decremented.
     */
    void setNode(int idx, MWNode<D, T> **node);

    /**
     * @brief Clear the node pointer in slot @p idx (set it to `nullptr`).
     * @param idx Linear box index.
     */
    void clearNode(int idx) { this->nodes[idx] = nullptr; }

    /**
     * @name Node access (mutable)
     * @{
     */

    /**
     * @brief Get the node stored at the box corresponding to @p idx.
     * @param idx Node index at the world scale.
     * @return Reference to the node.
     * @pre The slot must contain a non-null pointer.
     */
    MWNode<D, T> &getNode(NodeIndex<D> idx);

    /**
     * @brief Get the node stored at the box containing coordinate @p r.
     * @param r A point in world coordinates.
     * @return Reference to the node.
     * @pre The slot must contain a non-null pointer.
     */
    MWNode<D, T> &getNode(Coord<D> r);

    /**
     * @brief Get the node stored at linear index @p i.
     * @param i Linear box index (default 0).
     * @return Reference to the node.
     * @pre The slot must contain a non-null pointer.
     */
    MWNode<D, T> &getNode(int i = 0);
    ///@}

    /**
     * @name Node access (const)
     * @{
     */
    const MWNode<D, T> &getNode(NodeIndex<D> idx) const;
    const MWNode<D, T> &getNode(Coord<D> r) const;
    const MWNode<D, T> &getNode(int i = 0) const;
    ///@}

    /// @return Number of slots with non-null pointers.
    int getNOccupied() const { return this->nOccupied; }

    /// @return Raw pointer to the internal node-pointer array (size == size()).
    MWNode<D, T> **getNodes() { return this->nodes; }

protected:
    int nOccupied;        ///< Number of non-null entries in @ref nodes.
    MWNode<D, T> **nodes; ///< Dense array of node pointers (size equals number of boxes).

    /// Allocate and zero-initialize the @ref nodes array.
    void allocNodePointers();

    /// Clear all stored pointers (does not delete nodes).
    void deleteNodes();
};

} // namespace mrcpp
