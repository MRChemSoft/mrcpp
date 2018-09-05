#include "ProjectedNode.h"
#include "SerialTree.h"
#include "utils/Printer.h"

using namespace std;
using namespace Eigen;

namespace mrcpp {

template<int D>
void ProjectedNode<D>::createChildren() {
    MWNode<D>::createChildren();
    this->clearIsEndNode();
}

template<int D>
void ProjectedNode<D>::genChildren() {
    if (this->isBranchNode()) MSG_FATAL("Node already has children");
    this->tree->getSerialTree()->allocGenChildren(*this);
    this->setIsBranchNode();
}

template<int D>
void ProjectedNode<D>::deleteChildren() {
    MWNode<D>::deleteChildren();
    this->setIsEndNode();
}

template<int D>
void ProjectedNode<D>::dealloc() {
    int sIdx = this->serialIx;
    this->serialIx = -1;
    this->parentSerialIx = -1;
    this->childSerialIx = -1;
    this->tree->decrementNodeCount(this->getScale());
    this->tree->getSerialTree()->deallocNodes(sIdx);
}

/** Update the coefficients of the node by a mw transform of the scaling
  * coefficients of the children. Option to overwrite or add up existing
  * coefficients. Specialized for D=3 below. */
template<int D>
void ProjectedNode<D>::reCompress() {
    MWNode<D>::reCompress();
}

template<>
void ProjectedNode<3>::reCompress() {
    if (this->isBranchNode()) {
        if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");
	//can write directly from children coeff into parent coeff
	int stride = this->getMWChild(0).getNCoefs();
	double* inp  = this->getMWChild(0).getCoefs();
	double* out = this->coefs;
 
	assert(inp+7*stride == this->getMWChild(7).getCoefs());
 
	this->tree->getSerialTree()->S_mwTransformBack(inp, out, stride);
        this->setHasCoefs();
        this->calcNorms();
    }
}

template class ProjectedNode<1>;
template class ProjectedNode<2>;
template class ProjectedNode<3>;

} // namespace mrcpp
