/*
 *
 *
 *  \date Oct 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "GenNode.h"
#include "ProjectedNode.h"
#include "MWTree.h"

using namespace std;
using namespace Eigen;

template<int D>
GenNode<D>::GenNode(ProjectedNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx),
          genRootNode(&p) {
    this->allocCoefs(this->getTDim());
    this->zeroCoefs();
    this->setIsGenNode();
    this->clearHasWCoefs();
    this->tree->incrementGenNodeCount();
}

template<int D>
GenNode<D>::GenNode(GenNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx),
          genRootNode(p.genRootNode) {
    this->allocCoefs(this->getTDim());
    this->zeroCoefs();
    this->setIsGenNode();
    this->clearHasWCoefs();
    this->tree->incrementGenNodeCount();
}

template<int D>
GenNode<D>::~GenNode() {
    this->tree->decrementGenNodeCount(); //decrementNodeCount done in ~MWNode()
    if(this->tree->serialTree_p)this->tree->serialTree_p->DeAllocNodes(this->NodeRank);
    this->freeCoefs();
}

template<int D>
void GenNode<D>::createChild(int i) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::genChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MWNode<D> *child;
    if (this->tree->serialTree_p == 0){
      child = new GenNode<D>(*this, cIdx);
    } else {
      child = new (this->tree->serialTree_p->allocGenNodes(1))GenNode<D>(*this, cIdx);//GenNode also calls creator of MWNode
      child->NodeRank =  this->tree->serialTree_p->nNodes-1;
    }
    this->children[cIdx] = child;
}

template<int D>
void GenNode<D>::regenerateCoefs() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::allocCoefs(int nBlocks) {
    MWNode<D>::allocCoefs(nBlocks);
    this->tree->incrementAllocGenNodeCount();
}

template<int D>
void GenNode<D>::freeCoefs() {
    if (this->isAllocated()) {
        this->tree->decrementAllocGenNodeCount();
        MWNode<D>::freeCoefs();
    }
}

template<int D>
VectorXd& GenNode<D>::getCoefs() {
    lockSiblings();
    if (not this->hasCoefs()) {
        regenerateCoefs();
    }
    unlockSiblings();
    return MWNode<D>::getCoefs();
      //return new(&CoeffVector) Map<VectorXd>(*this->coefs,GenNodeNcoeff);

}

/** Get coefficients of GenNode, regenerate if needed, without locking. */
template<int D>
VectorXd& GenNode<D>::getCoefsNoLock() {
    if (not this->hasCoefs()) {
        regenerateCoefs();
    }
    return MWNode<D>::getCoefs();
}

template<int D>
const VectorXd& GenNode<D>::getCoefs() const {
    assert(this->hasCoefs());
    return MWNode<D>::getCoefs();
}

/** Clear coefficients of generated nodes.
  * The node structure is kept, only the coefficients are cleared. */
template<int D>
void GenNode<D>::clearGenerated() {
    this->freeCoefs();
    MWNode<D>::clearGenerated();
}

template<int D>
void GenNode<D>::lockSiblings() {
    MWNode<D> *parent = &this->getMWParent();
    if (parent != 0) {
	/* Since all threads set the locks in the same order starting from 0,
	there is no risk of a deadlock here. */
	for (int i = 0; i < parent->getNChildren(); i++) {
	    parent->getMWChild(i).lockNode();
	}
    }
}

template<int D>
void GenNode<D>::unlockSiblings() {
    MWNode<D> *parent = &this->getMWParent();
    if (parent != 0) {
	for (int i = 0; i < parent->getNChildren(); i++) {
	    parent->getMWChild(i).unlockNode();
	}
    }
}

template class GenNode<1> ;
template class GenNode<2> ;
template class GenNode<3> ;
