/**
 *
 *
 *  \date Aug 14, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 *
 */

#include "ProjectedNode.h"
#include "GenNode.h"
#include "TreeAllocator.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace std;
using namespace Eigen;

/** Root node constructor. By default root nodes are initialized to
  * represent functions which are constant zero. */
template<int D>
ProjectedNode<D>::ProjectedNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx)
        : FunctionNode<D> (t, nIdx) {
    this->allocCoefs(this->getTDim(), this->getKp1_d());
    this->setIsEndNode();
}

/** ProjectedNode constructor.
  * Creates an empty node given its parent and translation vector */
template<int D>
ProjectedNode<D>::ProjectedNode(ProjectedNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx) {
    this->allocCoefs(this->getTDim(), this->getKp1_d());
    this->setIsEndNode();
}

template<int D>
ProjectedNode<D>::~ProjectedNode() {
    this->freeCoefs();
}

/** Allocating child node.
  *
  * Given a child index, this routine creates an empty ProjectedNode child node
  * with the appropriate translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::createChild(int cIdx) {
    assert(this->children[cIdx] == 0);
//    ProjectedNode<D> *child;
//    if (this->tree->allocator == 0){
//        child = new ProjectedNode<D>(*this, cIdx);
//    } else {
//        child = new (this->tree->allocator->allocNodes(1))ProjectedNode<D>(*this, cIdx);
//        child = new ProjectedNode<D>(*this, cIdx);
//    }
    ProjectedNode<D> *child = new ProjectedNode<D>(*this, cIdx);
    this->children[cIdx] = child;
}

/** Generating child node.
  *
  * This routine creates 2^D children empty GenNodes with the appropriate
  * translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::genChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MWNode<D> *child = new GenNode<D>(*this, cIdx);
    this->children[cIdx] = child;
}

template class ProjectedNode<1> ;
template class ProjectedNode<2> ;
template class ProjectedNode<3> ;

