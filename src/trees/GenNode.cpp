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

#include "GenNode.h"
#include "NodeAllocator.h"
#include "utils/Printer.h"

namespace mrcpp {

template <int D> void GenNode<D>::createChildren() {
    NOT_REACHED_ABORT;
}

template <int D> void GenNode<D>::genChildren() {
    if (this->isBranchNode()) MSG_ABORT("Node already has children");
    this->getFuncTree().getGenNodeAllocator().allocChildren(*this);
    this->setIsBranchNode();
}

template <int D> void GenNode<D>::cvTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> void GenNode<D>::mwTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> void GenNode<D>::setValues(const Eigen::VectorXd &vec) {
    NOT_IMPLEMENTED_ABORT;
}

template <int D> void GenNode<D>::getValues(Eigen::VectorXd &vec) {
    MWNode<D> copy(*this);
    vec = Eigen::VectorXd::Zero(copy.getNCoefs());
    copy.mwTransform(Reconstruction);
    copy.cvTransform(Forward);
    for (int i = 0; i < this->n_coefs; i++) { vec(i) = copy.getCoefs()[i]; }
}

template <int D> double GenNode<D>::calcComponentNorm(int i) const {
    if (i == 0) {
        return MWNode<D>::calcComponentNorm(0);
    } else {
        return 0.0;
    }
}

template <int D> void GenNode<D>::dealloc() {
    int sIdx = this->serialIx;
    this->serialIx = -1;
    this->parentSerialIx = -1;
    this->childSerialIx = -1;
    this->tree->decrementGenNodeCount();
    this->getFuncTree().getGenNodeAllocator().deallocNodes(sIdx);
}

template <int D> void GenNode<D>::reCompress() {
    NOT_IMPLEMENTED_ABORT;
}

template class GenNode<1>;
template class GenNode<2>;
template class GenNode<3>;

} // namespace mrcpp
