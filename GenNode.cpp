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
    this->setIsGenNode();
    this->tree->incrementGenNodeCount();
}

template<int D>
GenNode<D>::GenNode(GenNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx),
          genRootNode(p.genRootNode) {
    this->setIsGenNode();
    this->tree->incrementGenNodeCount();
}

template<int D>
GenNode<D>::GenNode(const GenNode<D> &n)
        : FunctionNode<D> (n),
          genRootNode(n.genRootNode) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
GenNode<D>::~GenNode() {
    this->tree->decrementGenNodeCount();
    if (this->isAllocated()) {
        this->tree->decrementAllocGenNodeCount();
    }
}

template<int D>
void GenNode<D>::createChild(int i) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::genChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MWNode<D> *child = new GenNode<D>(*this, cIdx);
    this->children[cIdx] = child;
}

template<int D>
void GenNode<D>::regenerateCoefs() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::freeCoefs() {
    if (not this->isAllocated()) {
        return;
    }
    if (this->hasCoefs()) {
        this->calcNorms();
    }
    MWNode<D>::freeCoefs();
    this->tree->decrementAllocGenNodeCount();
}

template<int D>
void GenNode<D>::setCoefs(const VectorXd &c) {
    SET_NODE_LOCK();
    if (not this->isAllocated()) {
        this->getMWTree().incrementAllocGenNodeCount();
    }
    MWNode<D>::setCoefs(c);
    UNSET_NODE_LOCK();
}

template<int D>
VectorXd& GenNode<D>::getCoefs() {
    NOT_IMPLEMENTED_ABORT;
    lockSiblings();
    if (not this->hasCoefs()) {
        regenerateCoefs();
    }
    unlockSiblings();
    return MWNode<D>::getCoefs();
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
void GenNode<D>::mwTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::cvTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::lockSiblings() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void GenNode<D>::unlockSiblings() {
    NOT_IMPLEMENTED_ABORT;
}

template class GenNode<1> ;
template class GenNode<2> ;
template class GenNode<3> ;
