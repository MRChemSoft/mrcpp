/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "ProjectedNode.h"
#include "NodeAllocator.h"
#include "utils/Printer.h"
#include "utils/tree_utils.h"

using namespace Eigen;

namespace mrcpp {

template <int D> void ProjectedNode<D>::createChildren() {
    MWNode<D>::createChildren();
    this->clearIsEndNode();
}

template <int D> void ProjectedNode<D>::genChildren() {
    if (this->isBranchNode()) MSG_ABORT("Node already has children");
    this->getFuncTree().getGenNodeAllocator().allocChildren(*this);
    this->setIsBranchNode();
}

template <int D> void ProjectedNode<D>::deleteChildren() {
    MWNode<D>::deleteChildren();
    this->setIsEndNode();
}

template <int D> void ProjectedNode<D>::dealloc() {
    int sIdx = this->serialIx;
    this->serialIx = -1;
    this->parentSerialIx = -1;
    this->childSerialIx = -1;
    this->tree->decrementNodeCount(this->getScale());
    this->tree->getNodeAllocator().deallocNodes(sIdx);
}

/** Update the coefficients of the node by a mw transform of the scaling
 * coefficients of the children. Option to overwrite or add up existing
 * coefficients. Specialized for D=3 below. */
template <int D> void ProjectedNode<D>::reCompress() {
    MWNode<D>::reCompress();
}

template <> void ProjectedNode<3>::reCompress() {
    if (this->isBranchNode()) {
        if (not this->isAllocated()) MSG_ABORT("Coefs not allocated");
        // can write directly from children coeff into parent coeff
        int stride = this->getMWChild(0).getNCoefs();
        double *inp = this->getMWChild(0).getCoefs();
        double *out = this->coefs;

        assert(inp + 7 * stride == this->getMWChild(7).getCoefs());

        auto &tree = getMWTree();
        tree_utils::mw_transform_back(tree, inp, out, stride);
        this->setHasCoefs();
        this->calcNorms();
    }
}

template class ProjectedNode<1>;
template class ProjectedNode<2>;
template class ProjectedNode<3>;

} // namespace mrcpp
