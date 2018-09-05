#include "GenNode.h"
#include "SerialTree.h"
#include "utils/Printer.h"

namespace mrcpp {

template<int D>
void GenNode<D>::createChildren() {
    NOT_REACHED_ABORT;
}

template<int D>
void GenNode<D>::genChildren() {
    if (this->isBranchNode()) MSG_FATAL("Node already has children");
    this->tree->getSerialTree()->allocGenChildren(*this);
    this->setIsBranchNode();
}

template<int D>
void GenNode<D>::cvTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::mwTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::setValues(const Eigen::VectorXd &vec) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::getValues(Eigen::VectorXd &vec) {
    MWNode<D> copy(*this);
    vec = Eigen::VectorXd::Zero(copy.getNCoefs());
    copy.mwTransform(Reconstruction);
    copy.cvTransform(Forward);
    for (int i = 0; i < this->n_coefs; i++) {
        vec(i) = copy.getCoefs()[i];
    }
}

template<int D>
double GenNode<D>::calcComponentNorm(int i) const {
    if (i == 0) {
        return MWNode<D>::calcComponentNorm(0);
    } else {
        return 0.0;
    }
}

template<int D>
void GenNode<D>::dealloc() {
    int sIdx = this->serialIx;
    this->serialIx = -1;
    this->parentSerialIx = -1;
    this->childSerialIx = -1;
    this->tree->decrementGenNodeCount();
    this->tree->getSerialTree()->deallocGenNodes(sIdx);
}

template<int D>
void GenNode<D>::reCompress() {
    NOT_IMPLEMENTED_ABORT;
}

template class GenNode<1>;
template class GenNode<2>;
template class GenNode<3>;

} // namespace mrcpp
