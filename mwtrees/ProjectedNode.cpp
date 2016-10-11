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
#include "SerialTree.h"

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

    this->setIsLooseNode();//otherwise should not be allocated by constructor
    this->allocCoefs(this->getTDim(), this->getKp1_d());
    //this->tree->incrementNodeCount(this->getScale()); //loose node do not count

    this->setIsEndNode();
    this->setHasWCoefs();//default until known
}

/** ProjectedNode constructor.
  * Creates an empty node given its parent and translation vector */
template<int D>
ProjectedNode<D>::ProjectedNode(ProjectedNode<D> &p, int cIdx)
        : FunctionNode<D> (p, cIdx) {

    this->setIsLooseNode();//otherwise should not be allocated by constructor
    this->allocCoefs(this->getTDim(), this->getKp1_d());
    //this->tree->incrementNodeCount(this->getScale()); //loose node do not count

    this->setIsEndNode();
    this->setHasWCoefs();//default until known
}

template<int D>
ProjectedNode<D>::~ProjectedNode() {
  //decrementNodeCount done in ~MWNode
  if(this->tree->serialTree_p){
    if(this->isGenNode())cout<<"~ProjectedNode but Gen node flag "<<this->SNodeIx<<" "<<this->SNodeIx<<endl;
    if(not this->isLooseNode())cout<<"NOT LOOSE node! should not happen"<<endl;
    
    //DeAllocCoeff done in ~MWNode()
  }
}

/** Allocating child node.
  *
  * Given a child index, this routine creates an empty ProjectedNode child node
  * with the appropriate translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::createChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    ProjectedNode<D> *child;
    int NodeIx;
    double* coefs_p;
    if (this->tree->serialTree_p == 0){
        child = new ProjectedNode<D>(*this, cIdx);
	this->children[cIdx] = child;
    } else {
      MSG_FATAL("All ProjectedNodes siblings should be created together");
    }

}

/** Generating child node.
  *
  * This routine creates 2^D children empty GenNodes with the appropriate
  * translation and Hilbert path parameters. */
template<int D>
void ProjectedNode<D>::genChild(int cIdx) {
    assert(this->children[cIdx] == 0);
    MWNode<D> *child;
    int NodeIx;
    if (this->tree->serialTree_p == 0){
      child = new GenNode<D>(*this, cIdx);
    } else {
      MSG_FATAL("A ProjectedNodes should not be generated");
    }
    this->children[cIdx] = child;
}

template class ProjectedNode<1> ;
template class ProjectedNode<2> ;
template class ProjectedNode<3> ;

