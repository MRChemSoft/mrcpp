#include "TreeAllocator.h"
#include "MWNode.h"
#include "MWTree.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

/** FunctionTree constructor.
  * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
  * Initializes rootNodes to represent the zero function. */
template<int D>
TreeAllocator<D>::TreeAllocator(const MultiResolutionAnalysis<D> &mra,
                                int max_nodes)
        : nNodes(0),
          nNodesCoeff(0),
          lastNode(0),
          mwTree_p(0),
          dataArray(0),
          maxNodes(max_nodes),
          maxNodesCoeff(max_nodes),
          maxGenNodesCoeff(max_nodes/8),//to be thought of
          sizeTreeMeta(0),
          sizeNodeMeta(0),
          sizeNode(0) {

    this->sizeNodeMeta =8*((sizeof(ProjectedNode<D>)+7)/8);//we want only multiples of 8

    this->sizeNodeCoeff =(1<<D)*(MathUtils::ipow(mra.getOrder()+1,D))*sizeof(double) + sizeof(VectorXd);
    this->sizeGenNodeCoeff =(MathUtils::ipow(mra.getOrder()+1,D))*sizeof(double) + sizeof(VectorXd);
    println(0, "SizeNode Coeff (B) " << this->sizeNodeCoeff);
    println(0, "SizeNode Meta (B)  " << this->sizeNodeMeta*sizeof(double));

    //The first part of the Tree is filled with metadata; reserved size:
    this->sizeTreeMeta = 8*((sizeof(FunctionTree<D>)+7)/8);//we want only multiples of 8

    //The dynamical part of the tree is filled with nodes of size:
    this->sizeNode = this->sizeNodeMeta + sizeNodeCoeff;//in units of Bytes

    //Tree is defined as array of doubles, because C++ does not like void malloc
    this->dataArray = new double[this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNode/sizeof(double)];
    this->GenCoeffArray = new double[ this->maxGenNodesCoeff*this->sizeGenNodeCoeff/sizeof(double)];
    this->lastNode = (ProjectedNode<D>*) (this->dataArray + this->sizeTreeMeta/sizeof(double));
    this->lastNodeCoeff = (double*) (this->dataArray+this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNodeMeta/sizeof(double));

    //put a MWTree at start of array
    this->mwTree_p = new (this->dataArray) MWTree<D>(mra);

    this->mwTree_p->allocator = this;

    //reserve memory for root nodes
    //    ProjectedNode<D> *node_p = allocNodes(this->mwTree_p->getRootBox().size());

    //alloc root nodes
    NodeBox<D> &rBox = this->mwTree_p->getRootBox();
    MWNode<D> **roots = rBox.getNodes();
    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
        const NodeIndex<D> &nIdx = rBox.getNodeIndex(rIdx);
	//       roots[rIdx] = new (node_p) ProjectedNode<D>(*this->getTree(), nIdx);
        //node_p++;
         roots[rIdx] = new (allocNodes(1)) ProjectedNode<D>(*this->getTree(), nIdx);
    }
    this->mwTree_p->resetEndNodeTable();
}

//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* TreeAllocator<D>::allocNodes(int nAlloc) {
    this->nNodes += nAlloc;
    if (this->nNodes > this->maxNodes){
      println(0, "maxNodes exceeded " << this->maxNodes);
       this->nNodes -= nAlloc;
        return 0;
    } else {
      ProjectedNode<D>* oldlastNode  = this->lastNode;
      //We can have sizeNodeMeta with different size than ProjectedNode
      this->lastNode = (ProjectedNode<D>*) ((char *) this->lastNode + nAlloc*this->sizeNodeMeta);
        //println(0, "new size meta " << this->nNodes);
        return oldlastNode;
    }
}
//return pointer to the last active node or NULL if failed
template<int D>
GenNode<D>* TreeAllocator<D>::allocGenNodes(int nAlloc) {
    this->nNodes += nAlloc;
    if (this->nNodes > this->maxNodes){
      println(0, "maxNodes exceeded " << this->maxNodes);
       this->nNodes -= nAlloc;
        return 0;
    } else {
      GenNode<D>* oldlastNode  = (GenNode<D>*)this->lastNode;
      //We can have sizeNodeMeta with different size than ProjectedNode
      this->lastNode = (ProjectedNode<D>*) ((char *) this->lastNode + nAlloc*this->sizeNodeMeta);
        //println(0, "new size meta " << this->nNodes);
        return oldlastNode;
    }
}
//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* TreeAllocator<D>::allocCoeff(int nAllocCoeff) {
  //for now only scaing and wavelets, because we need to take into account the size of vector class
    if (nAllocCoeff!=8) MSG_FATAL("Only 2**D implemented now!");
    this->nNodesCoeff += 1;//nAllocCoeff;
    if (this->nNodesCoeff > this->maxNodesCoeff){
      println(0, "maxNodesCoeff exceeded " << this->maxNodesCoeff);
       this->nNodesCoeff -= 1;//nAllocCoeff;
        return 0;
    } else {
       double* oldlastNodeCoeff  = this->lastNodeCoeff;
      //We can have sizeNodeMeta with different size than ProjectedNode
       //       this->lastNodeCoeff = (double*) ((char *) this->lastNodeCoeff + nAllocCoeff*this->sizeNodeCoeff);
       this->lastNodeCoeff = (double*) ((char *) this->lastNodeCoeff + this->sizeNodeCoeff);
       //println(0, "new size Coeff " << this->nNodesCoeff);
      return oldlastNodeCoeff;
    }
}
//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* TreeAllocator<D>::allocGenCoeff(int nAllocCoeff) {
  //for now only scaing and wavelets, because we need to take into account the size of vector class
    if (nAllocCoeff!=8) MSG_FATAL("Only 2**D implemented now!");
    this->nGenNodesCoeff += 1;//nAllocCoeff;
    if (this->nGenNodesCoeff > this->maxGenNodesCoeff){
      println(0, "maxGenNodesCoeff exceeded " << this->maxGenNodesCoeff);
       this->nGenNodesCoeff -= 1;//nAllocCoeff;
        return 0;
    } else {
       double* oldlastGenNodeCoeff  = this->lastGenNodeCoeff;
      //We can have sizeNodeMeta with different size than ProjectedNode
       //       this->lastNodeCoeff = (double*) ((char *) this->lastNodeCoeff + nAllocCoeff*this->sizeNodeCoeff);
       this->lastGenNodeCoeff = (double*) ((char *) this->lastGenNodeCoeff + this->sizeGenNodeCoeff);
       //println(0, "new size Coeff " << this->nNodesCoeff);
      return oldlastGenNodeCoeff;
    }
}
/** TreeAllocator destructor. */
template<int D>
TreeAllocator<D>::~TreeAllocator() {
    MWNode<D> **roots = this->mwTree_p->getRootBox().getNodes();
    for (int i = 0; i < this->mwTree_p->getRootBox().size(); i++) {
        ProjectedNode<D> *node = static_cast<ProjectedNode<D> *>(roots[i]);
        node->~ProjectedNode();
        roots[i] = 0;
    }
    this->mwTree_p->~MWTree();
    delete[] this->dataArray;
    delete[] this->GenCoeffArray;
}

template class TreeAllocator<1>;
template class TreeAllocator<2>;
template class TreeAllocator<3>;
