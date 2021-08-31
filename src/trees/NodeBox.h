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

template <int D> class NodeBox final : public BoundingBox<D> {
public:
    NodeBox(const NodeIndex<D> &idx, const std::array<int, D> &nb = {});
    NodeBox(const NodeBox<D> &box);
    NodeBox(const BoundingBox<D> &box);
    NodeBox<D> &operator=(const NodeBox<D> &box) = delete;
    ~NodeBox() override;

    void setNode(int idx, MWNode<D> **node);
    void clearNode(int idx) { this->nodes[idx] = nullptr; }

    MWNode<D> &getNode(NodeIndex<D> idx);
    MWNode<D> &getNode(Coord<D> r);
    MWNode<D> &getNode(int i = 0);

    const MWNode<D> &getNode(NodeIndex<D> idx) const;
    const MWNode<D> &getNode(Coord<D> r) const;
    const MWNode<D> &getNode(int i = 0) const;

    int getNOccupied() const { return this->nOccupied; }
    MWNode<D> **getNodes() { return this->nodes; }

protected:
    int nOccupied;     ///< Number of non-zero pointers in box
    MWNode<D> **nodes; ///< Container of nodes

    void allocNodePointers();
    void deleteNodes();
};

} // namespace mrcpp
