#include "ProjectedNode.h"

using namespace std;
using namespace Eigen;

/** Update the coefficients of the node by a mw transform of the scaling
  * coefficients of the children. Option to overwrite or add up existing
  * coefficients. */
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
