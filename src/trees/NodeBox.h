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

#include "BoundingBox.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @class NodeBox
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 *
 * @brief Bounding box with node-pointer storage
 */
template <int D, typename T> class NodeBox final : public BoundingBox<D> {
public:
    /**
     * @brief Construct a NodeBox from a lower-corner index and number of boxes
     * @param idx Lower-corner @ref NodeIndex at the world scale
     * @param nb  Number of boxes per dimension (defaults to all ones)
     */
    NodeBox(const NodeIndex<D> &idx, const std::array<int, D> &nb = {});

    /**
     * @brief Copy-construct from another NodeBox
     * @param box Source NodeBox
     */
    NodeBox(const NodeBox<D, T> &box);

    /**
     * @brief Construct from a plain @ref BoundingBox
     * @param box Geometric box to take as base
     */
    NodeBox(const BoundingBox<D> &box);

    NodeBox<D, T> &operator=(const NodeBox<D, T> &box) = delete;

    /// @brief Destructor, deletes all nodes
    ~NodeBox() override;

    /**
     * @brief Store a node pointer in index @p idx
     * @param idx  Linear box index in `[0, size())`
     * @param node Address of the node pointer to store (double pointer)
     */
    void setNode(int idx, MWNode<D, T> **node);

    /**
     * @brief Clear the node pointer stored at index @p idx
     * @param idx Linear box index in `[0, size())`
     */
    void clearNode(int idx) { this->nodes[idx] = nullptr; }

    /**
     * @brief Get the node stored at the given index @p idx
     * @param idx Node index at the world scale
     * @return Reference to the node
     */
    MWNode<D, T> &getNode(NodeIndex<D> idx);
    /**
     * @brief Get the node stored at the box containing coordinate @p r.
     * @param r Coordinates of a point
     * @return Reference to the node
     */
    MWNode<D, T> &getNode(Coord<D> r);
    /**
     * @brief Get the node stored at the given index @p i
     * @param i Linear box index (default 0)
     * @return Reference to the node
     */
    MWNode<D, T> &getNode(int i = 0);

    /**
     * @brief Get the node stored at the given index @p idx
     * @param idx Node index at the world scale
     * @return Reference to the node
     */
    const MWNode<D, T> &getNode(NodeIndex<D> idx) const;
    /**
     * @brief Get the node stored at the box containing coordinate @p r.
     * @param r Coordinates of a point
     * @return Reference to the node
     */
    const MWNode<D, T> &getNode(Coord<D> r) const;
    /**
     * @brief Get the node stored at the given index @p i
     * @param i Linear box index (default 0)
     * @return Reference to the node
     */
    const MWNode<D, T> &getNode(int i = 0) const;

    int getNOccupied() const { return this->nOccupied; } ///< @return The number of occupied node slots
    MWNode<D, T> **getNodes() { return this->nodes; } ///< @return The nodes stored in this box

protected:
    int nOccupied;        ///< Number of non-null entries in @ref nodes.
    MWNode<D, T> **nodes; ///< Dense array of node pointers (size equals number of boxes).

    /// @brief Allocate the node double pointers
    void allocNodePointers();

    /// @brief Clear and delete all nodes
    void deleteNodes();
};

} // namespace mrcpp
