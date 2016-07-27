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
          maxNodesCoeff(max_nodes),//to be thought of
          maxGenNodesCoeff(max_nodes),//to be thought of
          sizeTreeMeta(0),
          sizeNodeMeta(0),
          sizeNode(0) {

    println(0, "max_nodes  " <<max_nodes);
    this->sizeNodeMeta =16*((sizeof(ProjectedNode<D>)+15)/16);//we want only multiples of 16

    this->sizeNodeCoeff =(1<<D)*(MathUtils::ipow(mra.getOrder()+1,D))*sizeof(double) + sizeof(VectorXd);
    //NB: Gen nodes should take less space?
    this->sizeGenNodeCoeff =(1<<D)*(MathUtils::ipow(mra.getOrder()+1,D))*sizeof(double) + sizeof(VectorXd);
    println(0, "SizeNode Coeff (B) " << this->sizeNodeCoeff);
    println(0, "SizeGenNode Coeff (B) " << this->sizeGenNodeCoeff);
    println(0, "SizeNode Meta (B)  " << this->sizeNodeMeta*sizeof(double));

    //The first part of the Tree is filled with metadata; reserved size:
    this->sizeTreeMeta = 16*((sizeof(FunctionTree<D>)+15)/16);//we want only multiples of 16

    //The dynamical part of the tree is filled with nodes of size:
    this->sizeNode = this->sizeNodeMeta + sizeNodeCoeff;//in units of Bytes

    //Tree is defined as array of doubles, because C++ does not like void malloc
    //NB: important to divide by sizeof(double) BEFORE multiplying. Otherwise we can get int overflow!
    int mysize = (this->sizeTreeMeta/sizeof(double) + this->maxNodes*(this->sizeNode/sizeof(double)));
    this->dataArray = new double[mysize];
    this->GenCoeffArray = new double[ this->maxGenNodesCoeff*this->sizeGenNodeCoeff/sizeof(double)];
    this->lastNode = (ProjectedNode<D>*) (this->dataArray + this->sizeTreeMeta/sizeof(double));
    //    this->lastNodeCoeff = (double*) (this->dataArray+this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNodeMeta/sizeof(double));
    this->lastNodeCoeff = (double*) (this->dataArray+this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNodeMeta/sizeof(double));//start after the metadata
    this->lastGenNodeCoeff = this->GenCoeffArray;//start at start of array
    this->nNodesCoeff=-1;//add 1 before each allocation
    this->nGenNodesCoeff=-1;//add 1 before each allocation

    println(0, "Allocated (MB)  " <<(mysize/1024/1024)*sizeof(double) );
    CoeffStack = new VectorXd * [maxNodesCoeff];
    CoeffStackStatus = new int[maxNodesCoeff+1];
    for (int i = 0; i <maxNodesCoeff;i++){
      CoeffStack[i]=new(this->lastNodeCoeff+i*this->sizeNodeCoeff/8) VectorXd((1<<D)*(MathUtils::ipow(mra.getOrder()+1,D)));
      CoeffStackStatus[i]=0;//0=unoccupied
    }
    CoeffStackStatus[maxNodesCoeff]=-1;//-1=unavailable
    //CoeffStackStatus[-1]=-1;//-1=unavailable TO FIX!

    GenCoeffStack = new VectorXd * [maxGenNodesCoeff];
    GenCoeffStackStatus = new int[maxGenNodesCoeff+1];
    for (int i = 0; i <maxGenNodesCoeff;i++){
      GenCoeffStack[i]=new(this->lastGenNodeCoeff+i*this->sizeGenNodeCoeff/sizeof(double)) VectorXd((1<<D)*(MathUtils::ipow(mra.getOrder()+1,D)));
      GenCoeffStackStatus[i]=0;//0=unoccupied
    }
    GenCoeffStackStatus[maxGenNodesCoeff]=-1;//-1=unavailable

     //put a MWTree at start of array
    this->mwTree_p = new (this->dataArray) MWTree<D>(mra);

    this->mwTree_p->allocator = this;

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
    //this->TempVector = new VectorXd((1<<D)*(MathUtils::ipow(mra.getOrder()+1,D)));
}

template<int D>
void TreeAllocator<D>::SerialTreeAdd(double c, FunctionTree<D>* &TreeB){
  println(0, " SerialTreeAdd ");

    int DepthMax = 100;
    int  slen = 0, counter = 0;
    int  slenA = 0, counterA = 0;
    MWNode<D>* stackA[DepthMax*8];
    int  slenB = 0, counterB = 0;
    MWNode<D>* stackB[DepthMax*8];
    FunctionTree<D>* TreeA=this->getTree();
    int N_GenCoeff=TreeA->getKp1_d();
    int N_Coeff=N_GenCoeff*TreeA->getTDim();
    //put root nodes on stack
    NodeBox<D> &rBoxA = TreeA->getRootBox();
    MWNode<D> **rootsA = rBoxA.getNodes();
    for (int rIdxA = 0; rIdxA < rBoxA.size(); rIdxA++) {
        const NodeIndex<D> &nIdx = rBoxA.getNodeIndex(rIdxA);
	stackA[slenA++] = TreeA->findNode(nIdx);
    }
    NodeBox<D> &rBoxB = TreeB->getRootBox();
    MWNode<D> **rootsB = rBoxB.getNodes();
    for (int rIdxB = 0; rIdxB < rBoxB.size(); rIdxB++) {
        const NodeIndex<D> &nIdx = rBoxB.getNodeIndex(rIdxB);
	stackB[slenB++] = TreeB->findNode(nIdx);
    }
 
    while (slenA) {
      counterA++;
      MWNode<D>* fposA = stackA[--slenA];
      // println(0,slenA<<" Children    " << fposA->getNChildren()<<" Node Depth  " <<fposA->getDepth());
	for (int i = 0; i < fposA->getNChildren(); i++) {
	  stackA[slenA++] = fposA->children[i];
	}
    }
    println(0, "TreeA Nodes   " <<counterA<<" Snodes "<<this->nNodes);

    while (slenB) {
      counterB++;
      MWNode<D>* fposB = stackB[--slenB];
      //            println(0,slenB<<" Children    " << fposB->getNChildren()<<" Node Depth  " <<fposB->getDepth());
	for (int i = 0; i < fposB->getNChildren(); i++) {
	  stackB[slenB++] = fposB->children[i];
	}
    }
    println(0, "TreeB Nodes   " <<counterB<<" Snodes "<<TreeB->allocator->nNodes);
    for (int rIdxA = 0; rIdxA < rBoxA.size(); rIdxA++) {
        const NodeIndex<D> &nIdx = rBoxA.getNodeIndex(rIdxA);
	stackA[slenA++] = TreeA->findNode(nIdx);
    }
    for (int rIdxB = 0; rIdxB < rBoxB.size(); rIdxB++) {
        const NodeIndex<D> &nIdx = rBoxB.getNodeIndex(rIdxB);
	stackB[slenB++] = TreeB->findNode(nIdx);
    }
 
    while (slenA+slenB) {
      counter++;
      MWNode<D>* fposA = stackA[--slenA];
      MWNode<D>* fposB = stackB[--slenB];
      //	println(0, slenA<<" TreeA children   " <<fposA->getNChildren()<<" TreeB children   " <<fposB->getNChildren());

      if(fposA->getNChildren()+fposB->getNChildren()){
	//println(0, slenA<<" TreeA children   " <<fposA->getNChildren()<<" TreeB children   " <<fposB->getNChildren());
	if(fposA->getNChildren()==0){
	  //println(0, slenA<<" AcreateChild   ");
	  for (int i = 0; i < fposA->getTDim(); i++)fposA->createChild(i);
	  fposA->setIsBranchNode();
	  //println(0, slenA<<" AgiveChildrenCoefs   ");
	  fposA->giveChildrenCoefs();
	//println(0, slenA<<" AgiveChildrenCoefs()succesful   ");
	}
	if(fposB->getNChildren()==0){
	  //println(0, slenB<<" BgenChildren()  ");
	  fposB->genChildren();
	  //println(0, slenB<<" BgenChildren() succesful ");
	}
	//	println(0, slenA<<" TreeA newchildren   " <<fposA->getNChildren()<<" TreeB newchildren   " <<fposB->getNChildren());

	for (int i = 0; i < fposA->getNChildren(); i++) {
	  stackA[slenA++] = fposA->children[i];
	  stackB[slenB++] = fposB->children[i];
	}
      }
      //      println(0, slenA<<" adding   "<<slenB);
      if(fposA->isGenNode()){
	for (int i=0; i<N_GenCoeff; i++)fposA->coefs[i]+=c*fposB->coefs[i];
	if(fposB->isGenNode()){
	  //should not happen	  
	  println(0, slenA<<" Adding two generated nodes? " );
	  //for (int i=N_GenCoeff; i<N_Coeff; i++)fposA->coefs[i]=0.0;
	}else{
	  for (int i=N_GenCoeff; i<N_Coeff; i++)fposA->coefs[i]=c*fposB->coefs[i];
	}	  
      }else{
	if(fposB->isGenNode()){
	  for (int i=0; i<N_GenCoeff; i++)fposA->coefs[i]+=c*fposB->coefs[i];
	}else{
	  for (int i=0; i<N_Coeff; i++)fposA->coefs[i]=c*fposB->coefs[i];
	}	  	
      }
      // println(0, slenA<<" added   "<<slenB);
    }
   println(0, "TreeAB Nodes   " <<counter);

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
      for (int i = 0; i<nAlloc; i++)oldlastNode->NodeRank=this->nNodes+i;
      return oldlastNode;
  }
}
template<int D>
void TreeAllocator<D>::DeAllocNodes(int NodeRank) {
  this->nNodes--;
  if (this->nNodes <0){
    println(0, "minNodes exceeded " << this->nNodes);
    this->nNodes++;
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
//double* TreeAllocator<D>::allocCoeff(int nAllocCoeff) {
VectorXd* TreeAllocator<D>::allocCoeff(int nAllocCoeff) {
  //for now only scaling and wavelets, because we need to take into account the size of vector class
  this->nNodesCoeff += 1;//nAllocCoeff;
    if (nAllocCoeff!=8) MSG_FATAL("Only 2**D implemented now!");
    if (this->nNodesCoeff > this->maxNodesCoeff ){
      println(0, "maxNodesCoeff exceeded " << this->maxNodesCoeff);
        return 0;
    } else if( CoeffStackStatus[this->nNodesCoeff]!=0){
      println(0, this->nNodesCoeff<<" CoeffStackStatus: not available " <<CoeffStackStatus[this->nNodesCoeff] );      
        return 0;
    }else{
      VectorXd* coeffVectorXd=this->CoeffStack[this->nNodesCoeff];
      this->CoeffStackStatus[this->nNodesCoeff]=1;
      //      println(0,this->nNodesCoeff<< " Coeff returned: " << coeffVectorXd);      
      return coeffVectorXd;
    }
}

//"Deallocate" the stack
template<int D>
void TreeAllocator<D>::DeAllocCoeff(int DeallocRank) {
    if (this->CoeffStackStatus[DeallocRank]==0)println(0, "deleting already unallocated coeff " << DeallocRank);
    this->CoeffStackStatus[DeallocRank]=0;//mark as available
    if(DeallocRank==this->nNodesCoeff){//top of stack
      int TopStack=this->nNodesCoeff;
      while(this->CoeffStackStatus[TopStack]==0 and TopStack>0){
	TopStack--;
      }
      this->nNodesCoeff=TopStack;//move top of stack
    }
}

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
VectorXd* TreeAllocator<D>::allocGenCoeff(int nAllocCoeff) {
  //for now only scaling and wavelets, because we need to take into account the size of vector class
    if (nAllocCoeff!=8) MSG_FATAL("Only 2**D implemented now!");
    this->nGenNodesCoeff += 1;//nAllocCoeff;
    if (this->nGenNodesCoeff > this->maxGenNodesCoeff ){
      println(0, "maxGenNodesCoeff exceeded " << this->maxGenNodesCoeff);
       this->nGenNodesCoeff -= 1;//nAllocCoeff;
        return 0;
    } else if( GenCoeffStackStatus[this->nGenNodesCoeff]!=0){
      println(0, this->nGenNodesCoeff<<" GenCoeffStackStatus: not available " <<GenCoeffStackStatus[this->nGenNodesCoeff] );      
        return 0;
    }else{
      VectorXd* coeffVectorXd=this->GenCoeffStack[this->nGenNodesCoeff];
      this->GenCoeffStackStatus[this->nGenNodesCoeff]=1;
      return coeffVectorXd;
    }
}

//"Deallocate" the stack
template<int D>
void TreeAllocator<D>::DeAllocGenCoeff(int DeallocRank) {
    if (this->GenCoeffStackStatus[DeallocRank]==0){
      println(0, "deleting already unallocated Gencoeff " << DeallocRank);
    }else{
      this->GenCoeffStackStatus[DeallocRank]=0;//mark as available
      if(DeallocRank==this->nGenNodesCoeff){//top of stack
	int TopStack=this->nGenNodesCoeff;
	while(this->CoeffStackStatus[TopStack]==0 and TopStack>0){
	  TopStack--;
	}
	this->nGenNodesCoeff=TopStack;//move top of stack
      }
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
    for (int i = 0; i <maxNodesCoeff;i++)this->CoeffStack[i]->~VectorXd();
    for (int i = 0; i <maxGenNodesCoeff;i++)this->GenCoeffStack[i]->~VectorXd();
    delete[] this->dataArray;
    delete[] this->GenCoeffArray;
    delete[] this->CoeffStack;
    delete[] this->CoeffStackStatus;
    delete[] this->GenCoeffStack;
    delete[] this->GenCoeffStackStatus;
}

template class TreeAllocator<1>;
template class TreeAllocator<2>;
template class TreeAllocator<3>;
