#include "TreeAllocator.h"
#include "MWNode.h"
#include "MWTree.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"

using namespace std;
using namespace Eigen;

/** FunctionTree constructor.
  * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
  * Initializes rootNodes to represent the zero function. */
template<int D>
TreeAllocator<D>::TreeAllocator(const MultiResolutionAnalysis<D> &mra, int maxNumberOfNodes)
        : nNodes(0),
          lastNode(0),
          mwTree_p(0),
          dataArray(0),
          maxNodes(maxNumberOfNodes),
          sizeTreeMeta(0),
          sizeNode(0) {
    NOT_IMPLEMENTED_ABORT;
    /*
    //The first part of the Tree is filled with metadata; reserved size:
    this->sizeTreeMeta = (sizeof(FunctionTree<D>)+7)/sizeof(double);
    //The dynamical part of the tree is filled with nodes of size:
    this->sizeNode = (sizeof(ProjectedNode<D>)+7)/sizeof(double);
    cout<<"SizeNode "<<sizeof(ProjectedNode<D>)<<endl;

    //Tree is defined as array of doubles, because C++ does not like void malloc
    this->dataArray = new double[this->sizeTreeMeta + maxNumberOfNodes*this->sizeNode];
    this->lastNode = (ProjectedNode<D>*) (this->dataArray + this->sizeTreeMeta);

    this->mwTree_p = new (this->dataArray) MWTree<D>(mra);//put a MWTree at start of array
    //this->funcTree_p = new FunctionTree<D>(mra);//put a MWTree at start of array

    ProjectedNode<D> *node_p = allocNodes(this->mwTree_p->getRootBox().size());

    NodeBox<D> &rBox = this->mwTree_p->getRootBox();
    MWNode<D> **roots = rBox.getNodes();
    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
        const NodeIndex<D> &nIdx = rBox.getNodeIndex(rIdx);
        roots[rIdx] = new (node_p) ProjectedNode<D>(getTree(), nIdx);
        node_p++;
    }
    this->mwTree_p->resetEndNodeTable();
//    this->mwTree_p->allocator = this;
    */
}

//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* TreeAllocator<D>::allocNodes(int nAlloc) {
    NOT_IMPLEMENTED_ABORT;
    /*
    this->nNodes += nAlloc;
    if (this->nNodes > this->maxNodes){
        this->nNodes -= nAlloc;
        return 0;
    } else {
        this->lastNode += nAlloc;
        cout << "new size " << this->nNodes << endl;
        return this->lastNode - nAlloc;
    }
    */
}
/** TreeAllocator destructor. */
template<int D>
TreeAllocator<D>::~TreeAllocator() {
    NOT_IMPLEMENTED_ABORT;
    /*
    MWNode<D> **roots = this->mwTree_p->getRootBox().getNodes();
    for (int i = 0; i < this->mwTree_p->getRootBox().size(); i++) {
        ProjectedNode<D> *node = static_cast<ProjectedNode<D> *>(roots[i]);
        node->~ProjectedNode();
        roots[i] = 0;
    }
    this->mwTree_p->~MWTree();
    delete[] this->dataArray;
    */
}

template class TreeAllocator<1> ;
template class TreeAllocator<2> ;
template class TreeAllocator<3> ;
