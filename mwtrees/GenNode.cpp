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
  
    this->setIsLooseNode();//otherwise should not be allocated by constructor
    this->allocCoefs(this->getTDim(), this->getKp1_d());
    this->zeroCoefs();
    this->setIsGenNode();
    this->clearHasWCoefs();
    this->tree->incrementGenNodeCount();
}

template<int D>
GenNode<D>::GenNode(GenNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx),
          genRootNode(p.genRootNode) {
    this->setIsLooseNode();//otherwise should not be allocated by constructor
    this->allocCoefs(this->getTDim(), this->getKp1_d());
    this->zeroCoefs();
    this->setIsGenNode();
    this->clearHasWCoefs();
    this->tree->incrementGenNodeCount();
}

template<int D>
GenNode<D>::~GenNode() {
    this->tree->decrementGenNodeCount(); //decrementNodeCount done in ~MWNode()
    //DeAllocGenCoeff done in ~MWNode()
    if(this->tree->serialTree_p){
      assert(this->isGenNode());
      assert(this->isLooseNode());
    }
}

template<int D>
void GenNode<D>::createChild(int i) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::genChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MWNode<D> *child;
    int NodeIx;
    if (this->tree->serialTree_p == 0){
      child = new GenNode<D>(*this, cIdx);
    } else {
    //NB: serial tree MUST generate all children consecutively
    //all children must be generated at once if several threads are active -> use genChildren
      MSG_FATAL("All GenNodes siblings in a Serial Tree Gen Nodes must be created at once");
      child = new (this->tree->serialTree_p->allocGenNodes(1, &NodeIx))GenNode<D>(*this, cIdx);
      child->NodeRank = NodeIx;
    }
    this->children[cIdx] = child;
}

template<int D>
void GenNode<D>::regenerateCoefs() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::allocCoefs(int n_blocks, int block_size) {
    MWNode<D>::allocCoefs(n_blocks, block_size);
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
void GenNode<D>::getCoefs(VectorXd &vec) {
    lockSiblings();
    if (not this->hasCoefs()) {
        regenerateCoefs();
    }
    unlockSiblings();

    MWNode<D>::getCoefs(vec);

}

/** Get coefficients of GenNode, regenerate if needed, without locking. */
template<int D>
void GenNode<D>::getCoefsNoLock(VectorXd &vec) {
    if (not this->hasCoefs()) {
        regenerateCoefs();
    }
    return MWNode<D>::getCoefs(vec);
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
