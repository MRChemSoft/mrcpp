#include <iostream>
#include "SerialTree.h"
#include "MWNode.h"
#include "MWTree.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"
#include "GenNode.h"
#include "MathUtils.h"
#include "Timer.h"
#include "parallel.h"

using namespace std;
using namespace Eigen;

/** SerialTree class constructor.
  * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
  * Initializes rootNodes to represent the zero function and allocate their nodes. 
  * NOTES:
  * Serial trees are made of projected nodes, and include gennodes and loose nodes separately.
  * All created (using class creator) Projected nodes or GenNodes are loose nodes. 
  * Loose nodes have their coeff in serial Tree, but not the node part. 
  * Projected nodes and GenNodes that are created by their creator, are detroyed by destructor ~ProjectedNode and ~GenNode. 
  * Serial tree nodes are not using the destructors, but explicitely call to deallocNodes or deallocGenNodes
  * Gen nodes and loose nodes are not counted with MWTree->[in/de]crementNodeCount()
*/
template<int D>
SerialTree<D>::SerialTree(FunctionTree<D> *tree, int max_nodes)
        : maxNodes(max_nodes),
          maxGenNodes(max_nodes),
          maxLooseNodesCoeff(64),
          nNodes(0),
          nGenNodes(0),
          nNodesCoeff(0),
          nGenNodesCoeff(0),
          nLooseNodesCoeff(0),
          tree_p(tree),
          lastNode(0),
          lastGenNode(0) {
    if (max_nodes == 0) this->maxLooseNodesCoeff = 0;

    //Size for GenNodes chunks. ProjectedNodes will be 8 times larger
    int sizePerChunk = 1024*1024;// 1 MB small for no waisting place, but large enough so that latency and overhead work is negligible

    this->sizeGenNodeCoeff = this->tree_p->getKp1_d();//One block
    this->sizeNodeCoeff = (1<<D)*this->sizeGenNodeCoeff;//TDim  blocks
    println(10, "SizeNode Coeff (kB) " << this->sizeNodeCoeff*sizeof(double)/1024);
    println(10, "SizeGenNode Coeff (kB) " << this->sizeGenNodeCoeff*sizeof(double)/1024);

    this->maxNodesPerChunk = sizePerChunk/this->sizeGenNodeCoeff;
    println(10, " max_nodes = " << max_nodes << ", nodes per chunk = " << this->maxNodesPerChunk);

    //indicate occupation of nodes
    this->nodeStackStatus = new int[this->maxNodes + 1];
    this->genNodeStackStatus = new int[this->maxGenNodes + 1];
    this->looseNodeStackStatus = new int[this->maxLooseNodesCoeff + 1];

    this->looseNodeCoeff = 0;//allocated when used first time
    this->looseCoeffStack = new double * [this->maxLooseNodesCoeff];


    this->lastNode = (ProjectedNode<D>*) this->sNodes;//position of last allocated node
    this->lastGenNode = this->sGenNodes;//position of last allocated Gen node

    this->nLooseNodesCoeff=0;

    //initialize stacks
    for (int i = 0; i < maxNodes; i++) {
        this->nodeStackStatus[i] = 0;//0=unoccupied
    }
    this->nodeStackStatus[maxNodes] = -1;//=unavailable

    //initialize stacks
    for (int i = 0; i < this->maxGenNodes; i++) {
        this->genNodeStackStatus[i] = 0;//0=unoccupied
    }
    this->genNodeStackStatus[this->maxGenNodes] = -1;//=unavailable

#ifdef HAVE_OPENMP
    omp_init_lock(&Stree_lock);
#endif

    //make virtual table pointers
    ProjectedNode<D>* tmpNode = new ProjectedNode<D>();
    this->cvptr_ProjectedNode =  *(char**)(tmpNode);
    delete tmpNode;

    GenNode<D>* tmpGenNode = new GenNode<D>();
    this->cvptr_GenNode =  *(char**)(tmpGenNode);
    delete tmpGenNode;

}

template<int D>
void SerialTree<D>::allocRoots(FunctionTree<D> &tree) {
    int sIx;
    double *coefs_p;
    //reserve place for nRoots
    int nRoots = tree.getRootBox().size();
    ProjectedNode<D> *root_p = this->allocNodes(nRoots, &sIx, &coefs_p);

    MWNode<D> **roots = tree.getRootBox().getNodes();
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
        roots[rIdx] = root_p;

        *(char**)(root_p) = this->cvptr_ProjectedNode;

        root_p->tree = &tree;
        root_p->parent = 0;
        for (int i = 0; i < (1 << D); i++) {
            root_p->children[i] = 0;
        }

        root_p->nodeIndex = tree.getRootBox().getNodeIndex(rIdx);
        root_p->hilbertPath = HilbertPath<D>();

        root_p->n_coefs = this->sizeNodeCoeff;
        root_p->coefs = coefs_p;

        root_p->serialIx = sIx;
        root_p->parentSerialIx = -1;//to indicate rootnode
        root_p->childSerialIx = -1;

        root_p->status = 0;

        root_p->clearNorms();
        root_p->setIsLeafNode();
        root_p->setIsAllocated();
        root_p->clearHasCoefs();
        root_p->setIsEndNode();
        root_p->setHasWCoefs();//default until known
        root_p->setIsRootNode();

        tree.incrementNodeCount(root_p->getScale());

#ifdef OPENMP
        omp_init_lock(&(root_p->node_lock));
#endif

        sIx++;
        root_p++;
        coefs_p += this->sizeNodeCoeff;
    }
}

template<int D>
void SerialTree<D>::allocChildren(FunctionNode<D> &parent) {
    int sIx;
    double *coefs_p;
    //NB: serial tree MUST generate all children consecutively
    //all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    ProjectedNode<D> *child_p = this->allocNodes(nChildren, &sIx, &coefs_p);

    //position of first child
    parent.childSerialIx = sIx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char**)(child_p) = this->cvptr_ProjectedNode;

        child_p->tree = parent.tree;
        child_p->parent = &parent;
        for (int i = 0; i < child_p->getTDim(); i++) {
            child_p->children[i] = 0;
        }

        child_p->nodeIndex = NodeIndex<D>(parent.getNodeIndex(), cIdx);
        child_p->hilbertPath = HilbertPath<D>(parent.getHilbertPath(), cIdx);

        child_p->n_coefs = this->sizeNodeCoeff;
        child_p->coefs = coefs_p;

        child_p->serialIx = sIx;
        child_p->parentSerialIx = parent.serialIx;
        child_p->childSerialIx = -1;

        child_p->status = 0;

        child_p->clearNorms();
        child_p->setIsLeafNode();
        child_p->setIsAllocated();
        child_p->clearHasCoefs();
        child_p->setIsEndNode();
        child_p->setHasWCoefs();

        child_p->tree->incrementNodeCount(child_p->getScale());

#ifdef OPENMP
        omp_init_lock(&child_p->node_lock);
#endif

        sIx++;
        child_p++;
        coefs_p += this->sizeNodeCoeff;
    }
}

template<int D>
void SerialTree<D>::allocGenChildren(FunctionNode<D> &parent) {
    NOT_IMPLEMENTED_ABORT;

    /*
    if (this->isBranchNode()) MSG_FATAL("Node already has children");
    int nChildren = this->getTDim();

    //NB: serial tree MUST generate all children consecutively
    //all children must be generated at once if several threads are active
    if (this->tree->serialTree_p){
      int NodeIx;
      double *coefs_p;
      //reserve place for nChildren
      GenNode<D>* GenNode_p = this->tree->serialTree_p->allocGenNodes(nChildren, &NodeIx, &coefs_p);
      MWNode<D> *child;
      this->childSerialIx = NodeIx;//not used fro Gennodes?
      for (int cIdx = 0; cIdx < nChildren; cIdx++) {
	//if(true){
	if(false){
          NOT_IMPLEMENTED_ABORT;
	  if (ProjectedNode<D> *node = dynamic_cast<ProjectedNode<D> *>(this)) {
	    child = new (GenNode_p)GenNode<D>(*node, cIdx);
	  } else if (GenNode<D> *node = dynamic_cast<GenNode<D> *>(this)){
	    child = new (GenNode_p)GenNode<D>(*node, cIdx);
	  }else{
	    MSG_FATAL("genChildren error");
	  }
	  this->children[cIdx] = child;
	  //Noderank is written in allocGenNodes already!
	  //child->serialIx = NodeIx;
	}else{

 	  *(char**)(GenNode_p)=this->tree->serialTree_p->cvptr_GenNode;
	  GenNode_p->serialIx = NodeIx;
	  GenNode_p->tree = this->tree;
	  GenNode_p->parent = this;
	  GenNode_p->parentSerialIx = this->serialIx;
	  GenNode_p->nodeIndex = NodeIndex<D>(this->getNodeIndex(),cIdx);
	  GenNode_p->hilbertPath = HilbertPath<D>(this->getHilbertPath(), cIdx);
	  GenNode_p->squareNorm = -1.0;
	  GenNode_p->status = 0;

	  GenNode_p->zeroNorms();
	  GenNode_p->childSerialIx = -1;
	  for (int i = 0; i < GenNode_p->getTDim(); i++) {
	    GenNode_p->children[i] = 0;
	  }
	  GenNode_p->setIsLeafNode();
	  GenNode_p->coefs = coefs_p;
	  GenNode_p->n_coefs = this->tree->serialTree_p->sizeGenNodeCoeff;
	  GenNode_p->setIsAllocated();
	  GenNode_p->setHasCoefs();

	  //	  GenNode_p->zeroCoefs();//NB: not initialized

	  GenNode_p->tree->incrementGenNodeCount();
	  GenNode_p->setIsGenNode();
	  
	  this->children[cIdx] = GenNode_p;

#ifdef OPENMP
	  omp_init_lock(&GenNode_p->node_lock);
#endif

	}

	GenNode_p++;
	NodeIx++;
        coefs_p += this->tree->serialTree_p->sizeGenNodeCoeff;	
	
      }
    }else{
      NOT_IMPLEMENTED_ABORT;
      //for (int cIdx = 0; cIdx < nChildren; cIdx++) {
	//genChild(cIdx);
      //}
    }
    this->setIsBranchNode();
    */
}

/** Overwrite all pointers defined in the tree.
  * Necessary after sending the tree 
  * could be optimized. Should reset other counters? (GenNodes...) */
template<int D>
void SerialTree<D>::rewritePointers(int nChunks){
    NOT_IMPLEMENTED_ABORT;
    /*
  int depthMax = 100;
  MWNode<D>* stack[depthMax*8];
  int slen = 0, counter = 0;

  this->nGenNodes = 0;
  this->nGenNodesCoeff = -1;
  this->nLooseNodesCoeff = 0;

  //reinitialize stacks
  for (int i = 0; i < this->maxNodes; i++) {
        this->nodeStackStatus[i] = 0;
  }

  for (int i = 0; i < this->maxGenNodes; i++) {
        this->genNodeStackStatus[i] = 0;//0=unoccupied
  }
  this->genNodeStackStatus[this->maxGenNodes] = -1;//=unavailable

  for (int i = 0; i < this->maxLooseNodesCoeff; i++) {
    this->looseCoeffStackStatus[i] = 0;//0=unoccupied
  }
  this->looseCoeffStackStatus[this->maxLooseNodesCoeff]=-1;//-1=unavailable

  this->getTree()->nNodes = 0;
  this->getTree()->nodesAtDepth.clear();
  this->getTree()->squareNorm = 0.0;

  for(int ichunk = 0 ; ichunk < nChunks; ichunk++){
    for(int inode = 0 ; inode < this->maxNodesPerChunk; inode++){
      ProjectedNode<D>* node = (this->nodeChunks[ichunk]) + inode;
      if (node->SNodeIx >= 0) {
	  //node is part of tree, should be processed
	  this->getTree()->incrementNodeCount(node->getScale());
	  if (node->isEndNode()) this->getTree()->squareNorm += node->getSquareNorm();
	  
	  //normally (intel) the virtual table does not change, but we overwrite anyway
	  *(char**)(Node) = this->cvptr_ProjectedNode;
	  
	  Node->tree = this->getTree();

	  //"adress" of coefs is the same as node, but in another array
	  Node->coefs = this->NodeCoeffChunks[ichunk]+ inode*this->sizeNodeCoeff;
	  
	  //adress of parent and children must be corrected
	  //can be on a different chunks
	  if(Node->parentSNodeIx>=0){
	    int n_ichunk = Node->parentSNodeIx/this->maxNodesPerChunk;
	    int n_inode = Node->parentSNodeIx%this->maxNodesPerChunk;
	    Node->parent = this->NodeChunks[n_ichunk] + n_inode;
	  }else{Node->parent = 0;}
	    
	  for (int i = 0; i < Node->getNChildren(); i++) {
	    int n_ichunk = (Node->childSNodeIx+i)/this->maxNodesPerChunk;
	    int n_inode = (Node->childSNodeIx+i)%this->maxNodesPerChunk;
	    Node->children[i] = this->NodeChunks[n_ichunk] + n_inode;
	  }
	  this->NodeStackStatus[Node->SNodeIx] = 1;//occupied
#ifdef OPENMP
	  omp_init_lock(&(Node->node_lock));
#endif
	}

    }
  }

  //update other MWTree data
  FunctionTree<D>* Tree = static_cast<FunctionTree<D>*> (this->mwTree_p);

  NodeBox<D> &rBox = Tree->getRootBox();
  MWNode<D> **roots = rBox.getNodes();

  for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
    roots[rIdx] = (this->NodeChunks[0]) + rIdx;//adress of roots are at start of NodeChunks[0] array
  }

  this->getTree()->resetEndNodeTable();

*/
}

/** Make 8 children nodes with scaling coefficients from parent
 * Does not put 0 on wavelets
 */
template<int D>
void SerialTree<D>::GenS_nodes(MWNode<D>* Node){
    NOT_IMPLEMENTED_ABORT;
/*
  bool ReadOnlyScalingCoeff=true;

  //  if(Node->isGenNode())ReadOnlyScalingCoeff=true;
  if(Node->hasWCoefs())ReadOnlyScalingCoeff=false;
  double* cA;

  Node->genChildren();//will make children and allocate coeffs, but without setting values for coeffs.


  double* coeffin  = Node->coefs;
  double* coeffout = Node->children[0]->coefs;

  int Children_Stride = this->sizeGenNodeCoeff;
  S_mwTransform(coeffin, coeffout, ReadOnlyScalingCoeff, Children_Stride);
*/
}


/** Make children scaling coefficients from parent
 * Other node info are not used/set
 * coeff_in are not modified.
 * The output is written directly into the 8 children scaling coefficients. 
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
template<int D>
void SerialTree<D>::S_mwTransform(double* coeff_in, double* coeff_out, bool ReadOnlyScalingCoeff, int Children_Stride, bool b_overwrite) {
    NOT_IMPLEMENTED_ABORT;
    /*

  int operation = Reconstruction;
  int kp1 = this->getTree()->getKp1();
  int tDim = (1<<D);
  int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
  int kp1_d = kp1_dm1*kp1;//
  const MWFilter &filter = this->getTree()->getMRA().getFilter();
  double overwrite = 0.0;
  double *tmp;
  double tmpcoeff[kp1_d*tDim];
  double tmpcoeff2[kp1_d*tDim];
  int ftlim=tDim;
  int ftlim2=tDim;
  int ftlim3=tDim;
  if(ReadOnlyScalingCoeff){
    ftlim=1;
    ftlim2=2;
    ftlim3=4;
    //NB: Careful: tmpcoeff tmpcoeff2 are not initialized to zero
    //must not read these unitialized values!
  }

  overwrite = 0.0;
  int i = 0;
  int mask = 1;
  for (int gt = 0; gt < tDim; gt++) {
    double *out = tmpcoeff + gt * kp1_d;
    for (int ft = 0; ft < ftlim; ft++) {
      // Operate in direction i only if the bits along other
      // directions are identical. The bit of the direction we
      // operate on determines the appropriate filter/operator
      if ((gt | mask) == (ft | mask)) {
	double *in = coeff_in + ft * kp1_d;
	int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	//overwrite = false;
	overwrite = 1.0;
      }
    }
    //overwrite = true;
    overwrite = 0.0;
  }
  if(D>1){
    i++;
    mask = 2;//1 << i;
    for (int gt = 0; gt < tDim; gt++) {
      double *out = tmpcoeff2 + gt * kp1_d;
      for (int ft = 0; ft < ftlim2; ft++) {
      // Operate in direction i only if the bits along other
      // directions are identical. The bit of the direction we
      // operate on determines the appropriate filter/operator
	if ((gt | mask) == (ft | mask)) {
	  double *in = tmpcoeff + ft * kp1_d;
	  int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	  const MatrixXd &oper = filter.getSubFilter(filter_index, operation);
	  
	MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	overwrite = 1.0;
	}
      }
      overwrite = 0.0;
    }
  }
  if(D>2){
    overwrite = 1.0;
    if(b_overwrite)overwrite = 0.0;
    i++;
    mask = 4;//1 << i;
    for (int gt = 0; gt < tDim; gt++) {
      //double *out = coeff_out + gt * kp1_d;
      double *out = coeff_out + gt * Children_Stride;//write right into children
      for (int ft = 0; ft < ftlim3; ft++) {
	// Operate in direction i only if the bits along other
	// directions are identical. The bit of the direction we
	// operate on determines the appropriate filter/operator
      if ((gt | mask) == (ft | mask)) {
	double *in = tmpcoeff2 + ft * kp1_d;
	int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	const MatrixXd &oper = filter.getSubFilter(filter_index, operation);
	
	MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	overwrite = 1.0;
      }
      }
      overwrite = 1.0;
      if(b_overwrite)overwrite = 0.0;
    }
  }

  if(D>3)MSG_FATAL("D>3 NOT IMPLEMENTED for S_mwtransform");

  if(D<3){
    double *out;
    if(D==1)out=tmpcoeff;
    if(D==2)out=tmpcoeff2;
    if(b_overwrite){
      for (int j = 0; j < tDim; j++){ 
	for (int i = 0; i < kp1_d; i++){ 
	  coeff_out[i+j*Children_Stride]=out[i+j*kp1_d];
	}
      }
    }else{
      for (int j = 0; j < tDim; j++){ 
	for (int i = 0; i < kp1_d; i++){ 
	  coeff_out[i+j*Children_Stride]+=out[i+j*kp1_d];
	}
      }
    }
  }
  */
}

/** Make parent from children scaling coefficients
 * Other node info are not used/set
 * coeff_in are not modified.
 * The output is read directly from the 8 children scaling coefficients. 
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
template<int D>
void SerialTree<D>::S_mwTransformBack(double* coeff_in, double* coeff_out, int Children_Stride) {
  assert(D==3);
  int operation = Compression;
  int kp1 = this->getTree()->getKp1();
  int kp1_d = this->getTree()->getKp1_d();
  int tDim = (1<<D);
  int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
  const MWFilter &filter = this->getTree()->getMRA().getFilter();
  double overwrite = 0.0;
  double tmpcoeff[kp1_d*tDim];

  int ftlim = tDim;
  int ftlim2 = tDim;
  int ftlim3 = tDim;

  int i = 0;
  int mask = 1;
  for (int gt = 0; gt < tDim; gt++) {
        double *out = coeff_out + gt * kp1_d;
        for (int ft = 0; ft < ftlim; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
	        double *in = coeff_in + ft * Children_Stride;
	        int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	        const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	        MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	        overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    i++;
    mask = 2;//1 << i;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = tmpcoeff + gt * kp1_d;
        for (int ft = 0; ft < ftlim2; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
	        double *in = coeff_out + ft * kp1_d;
	        int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	        const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	        MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	        overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    i++;
    mask = 4;//1 << i;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = coeff_out + gt * kp1_d;
        //double *out = coeff_out + gt * N_coeff;
        for (int ft = 0; ft < ftlim3; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
	        double *in = tmpcoeff + ft * kp1_d;
	        int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	        const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	        MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	        overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
  //  if(D%2)for(int i=0; i<kp1_d*tDim; i++) coeff_out[i] = tmpcoeff[i];
}


/** Allocating a Projected Serial (root) node.
  *
  * This routine creates an empty ProjectedNode node
  * with the appropriate translation */
template<int D>
ProjectedNode<D>* SerialTree<D>::createSnode(const NodeIndex<D> & nIdx) {
    NOT_IMPLEMENTED_ABORT;
/*

  int NodeIx;
  double *coefs_p;
  ProjectedNode<D>* newNode=this->allocNodes(1, &NodeIx, &coefs_p);

  *(char**)(newNode) = this->cvptr_ProjectedNode;
  newNode->SNodeIx = NodeIx;
  newNode->tree = this->mwTree_p;
  newNode->parent = 0;
  newNode->parentSNodeIx = -1;//to indicate rootnode
  newNode->nodeIndex = nIdx;
  newNode->hilbertPath = HilbertPath<D>();
  newNode->squareNorm = -1.0;
  newNode->status = 0;
  newNode->clearNorms();
  newNode->childSNodeIx = -1;
  for (int i = 0; i < (1 << D); i++) {
    newNode->children[i] = 0;
  }
  newNode->setIsLeafNode();
  newNode->coefs = coefs_p;
  newNode->n_coefs = this->sizeNodeCoeff;
  newNode->setIsAllocated();
  newNode->clearHasCoefs();
  
  newNode->tree->incrementNodeCount(newNode->getScale());
  newNode->setIsEndNode();
  newNode->setHasWCoefs();//default until known
  
#ifdef OPENMP
  omp_init_lock(&(newNode->node_lock));
#endif

  return newNode;
*/
}


//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* SerialTree<D>::allocNodes(int nAlloc, int* serialIx, double** coefs_p) {
    *serialIx = this->nNodes;
    int chunkIx = *serialIx%(this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx+nAlloc > this->maxNodesPerChunk ) {
        if (this->nNodes+nAlloc >= this->maxNodes){
            println(0, "maxNodes exceeded " << this->maxNodes);
            MSG_FATAL("maxNodes exceeded ");
        }

        //we want nodes allocated simultaneously to be allocated on the same pice.
        //possibly jump over the last nodes from the old chunk
        this->nNodes = this->maxNodesPerChunk*((this->nNodes+nAlloc-1)/this->maxNodesPerChunk);//start of next chunk

        int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk

        //careful: NodeChunks.size() is an unsigned int
        if (chunk+1 > this->nodeChunks.size()){
	    //need to allocate new chunk
	    this->sNodes = (ProjectedNode<D>*) new char[this->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
	    this->nodeChunks.push_back(this->sNodes);
            this->sNodesCoeff = new double[this->sizeNodeCoeff*this->maxNodesPerChunk];
            this->nodeCoeffChunks.push_back(this->sNodesCoeff);
	    if (chunk%100==99 and D==3) println(10,endl<<" number of nodes "<<this->nNodes <<",number of Nodechunks now " << this->nodeChunks.size()<<", total size coeff  (MB) "<<(this->nNodes/1024) * this->sizeNodeCoeff/128);
        }
        this->lastNode = this->nodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);
        *serialIx = this->nNodes;
        chunkIx = *serialIx%(this->maxNodesPerChunk);
    }
    assert((this->nNodes+nAlloc-1)/this->maxNodesPerChunk < this->nodeChunks.size());

    ProjectedNode<D> *newNode  = this->lastNode;
    ProjectedNode<D> *newNode_cp  = newNode;
    *coefs_p = this->sNodesCoeff + chunkIx*this->sizeNodeCoeff;
 
    //We can have sizeNodeMeta with different size than ProjectedNode
    //println(0, "new size meta " << this->nNodes);
    for (int i = 0; i < nAlloc; i++) {
        if (this->nodeStackStatus[*serialIx+i] != 0)
	    println(0, *serialIx+i<<" NodeStackStatus: not available " << this->nodeStackStatus[*serialIx+i]);
        this->nodeStackStatus[*serialIx+i] = 1;
        newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->lastNode += nAlloc;
    return newNode;
}

template<int D>
void SerialTree<D>::deallocNodes(int serialIx) {
    if (this->nNodes < 0) {
        println(0, "minNodes exceeded " << this->nNodes);
        this->nNodes++;
    }
    this->nodeStackStatus[serialIx] = 0;//mark as available
    if(serialIx == this->nNodes-1){//top of stack
        int topStack = this->nNodes;
        while (this->nodeStackStatus[topStack-1] == 0){
            topStack--;
            if (topStack < 1) break;
        }
        this->nNodes = topStack;//move top of stack
        //has to redefine lastGenNode
        int chunk = this->nNodes/this->maxNodesPerChunk;//find the right chunk
        this->lastNode = this->nodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);
    }
}

//return pointer to the last active node or NULL if failed
template<int D>
GenNode<D>* SerialTree<D>::allocGenNodes(int nAlloc, int* NodeIx, double** coefs_p) {
    NOT_IMPLEMENTED_ABORT;
/*
  omp_set_lock(&Stree_lock);
  *NodeIx = this->nGenNodes;
  int ChunkIx = *NodeIx%(this->maxNodesPerChunk);
  
  //Not necessarily wrong, but new:
  assert(nAlloc == (1<<D));

  if(ChunkIx == 0  or  ChunkIx+nAlloc > maxNodesPerChunk ){
    //start on new chunk

      if (this->nGenNodes+nAlloc >= this->maxGenNodes){
	println(0, "maxNodes exceeded " << this->maxGenNodes);
	MSG_FATAL("maxNodes exceeded ");
      } 

      //we want nodes allocated simultaneously to be allocated on the same chunk.
      //possibly jump over the last nodes from the old chunk
      this->nGenNodes=this->maxNodesPerChunk*((this->nGenNodes+nAlloc-1)/maxNodesPerChunk);//start of next chunk
     //
      int chunk = this->nGenNodes/maxNodesPerChunk;//find the right chunk
    //careful: NodeChunks.size() is an unsigned int
      if(chunk+1 > this->GenNodeChunks.size()){
	//need to allocate new chunk
	this->SGenNodes = (GenNode<D>*) new char[this->maxNodesPerChunk*sizeof(GenNode<D>)];
	this->GenNodeChunks.push_back(this->SGenNodes);
	this->SGenNodesCoeff = new double[this->sizeGenNodeCoeff*this->maxNodesPerChunk];
	
	GenNodeCoeffChunks.push_back(this->SGenNodesCoeff);
	if(chunk%100==99 and D==3)println(10,endl<<" number of Gennodes "<<this->nGenNodes <<",number of GenNodechunks now " << this->GenNodeChunks.size()<<", total size coeff  (MB) "<<(this->nGenNodes/1024) * this->sizeGenNodeCoeff/128);

      }
      this->lastGenNode = this->GenNodeChunks[chunk] + this->nGenNodes%(this->maxNodesPerChunk);
      
      *NodeIx = this->nGenNodes; 
      ChunkIx = *NodeIx%(this->maxNodesPerChunk);
      
    }
    
    assert((this->nGenNodes+nAlloc-1)/maxNodesPerChunk < this->GenNodeChunks.size());

    GenNode<D>* newNode  = this->lastGenNode;
    GenNode<D>* newNode_cp  = newNode;
    *coefs_p = this->SGenNodesCoeff + ChunkIx*this->sizeGenNodeCoeff;
    for (int i = 0; i<nAlloc; i++){
      newNode_cp->SNodeIx = *NodeIx+i;//Until overwritten!
      if (this->GenNodeStackStatus[*NodeIx+i]!=0)
	println(0, *NodeIx+i<<" NodeStackStatus: not available " << this->GenNodeStackStatus[*NodeIx+i]);
      this->GenNodeStackStatus[*NodeIx+i]=1;
      newNode_cp++;
    }
    
    this->nGenNodes += nAlloc;
    this->lastGenNode += nAlloc;
    
    omp_unset_lock(&Stree_lock);
    return newNode;
    */
}

template<int D>
void SerialTree<D>::deallocGenNodes(int SNodeIx) {
    NOT_IMPLEMENTED_ABORT;
/*
  omp_set_lock(&Stree_lock);
  if (this->nGenNodes <0){
    println(0, "minNodes exceeded " << this->nGenNodes);
    this->nGenNodes++;
  }
  this->GenNodeStackStatus[SNodeIx]=0;//mark as available
  if(SNodeIx==this->nGenNodes-1){//top of stack
    int TopStack=this->nGenNodes;
    while(this->GenNodeStackStatus[TopStack-1]==0){
      TopStack--;
      if(TopStack<1)break;
    }

    this->nGenNodes=TopStack;//move top of stack
    //has to redefine lastGenNode
    int chunk = this->nGenNodes/maxNodesPerChunk;//find the right chunk
    //    cout<<chunk<<" "<<this->GenNodeChunks.size()<<" old lastgennode "<<(double*) this->lastGenNode;
    this->lastGenNode = this->GenNodeChunks[chunk] + this->nGenNodes%(this->maxNodesPerChunk);
    //    cout<<" new lastgennode "<<(double*) this->lastGenNode<<endl;

  }

  omp_unset_lock(&Stree_lock);
  */
 }

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocCoeff(int nAllocCoeff, MWNode<D>* node) {
    NOT_IMPLEMENTED_ABORT;
/*
  if(node->isLooseNode()){
    return this->allocLooseCoeff(nAllocCoeff, node);
  }else{
     MSG_FATAL("ERROR SHOULD NOT ALLOCATE COEFF HERE " );
  }
  */
}

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocLooseCoeff(int nAllocCoeff, MWNode<D>* node) {
    NOT_IMPLEMENTED_ABORT;
/*

  if (nAllocCoeff!=1<<D) MSG_FATAL("Only 2**D implemented now!");

  //Each omp thread use own part of array, no locks!

  if( this->LooseNodeCoeff==0){
    if(omp_get_thread_num()==0){
      double * LooseNodesCoeff_p = new double[this->maxLooseNodesCoeff*this->sizeNodeCoeff];
      for (int i = 0; i <this->maxLooseNodesCoeff;i++){
	this->LooseCoeffStack[i] = LooseNodesCoeff_p+i*this->sizeNodeCoeff;
	this->LooseCoeffStackStatus[i] = 0;//0=unoccupied
      }
      this->LooseCoeffStackStatus[this->maxLooseNodesCoeff]=-1;//-1=unavailable
      this->LooseNodeCoeff = LooseNodesCoeff_p;//MUST be set after Stack inits
    }else{
      //wait until allocated
      while(this->LooseNodeCoeff==0){usleep(1000);}
    }
  }

  if (this->nLooseNodesCoeff+1 >= this->maxLooseNodesCoeff ){
    println(0, "maxLooseNodesCoeff exceeded " << this->maxLooseNodesCoeff);
    MSG_FATAL("maxLooseNodesCoeff exceeded ");
       // omp_unset_lock(&Stree_lock);
    return 0;
  } else if( LooseCoeffStackStatus[this->nLooseNodesCoeff]!=0 and omp_get_num_threads()==1){
    println(0, this->nLooseNodesCoeff<<" LooseCoeffStackStatus: not available " <<LooseCoeffStackStatus[this->nLooseNodesCoeff] <<" tree at "<<this->mwTree_p);     
    return 0;
  }else{
    int myindex = this->nLooseNodesCoeff;//default for one thread
    if(omp_get_num_threads()>1){
      //find the lowest available index that are on own part of stack
      myindex=omp_get_thread_num();
      while(this->LooseCoeffStackStatus[myindex]==1){
	myindex += omp_get_num_threads();
	if(myindex>this->maxLooseNodesCoeff)MSG_FATAL("maxLooseNodesCoeff exceeded for one thread ");	
      }
   
    }
    double* LooseNodeCoeff=this->LooseCoeffStack[myindex];
    this->LooseCoeffStackStatus[myindex]=1;
    node->SNodeIx = myindex;
#pragma omp atomic
    (this->nLooseNodesCoeff)++;

    return LooseNodeCoeff;
  }
*/
}


//"Deallocate" the stack
template<int D>
void SerialTree<D>::deallocLooseCoeff(int DeallocIx) {
    NOT_IMPLEMENTED_ABORT;
/*
  if (this->LooseCoeffStackStatus[DeallocIx]==0){
    println(0, "deleting already unallocated loose coeff " << DeallocIx);
  }
  this->LooseCoeffStackStatus[DeallocIx]=0;//mark as available
  if(omp_get_num_threads()==1){
  if(DeallocIx==this->nLooseNodesCoeff-1){//top of stack
    int TopStack=this->nLooseNodesCoeff-1;
    while(this->LooseCoeffStackStatus[TopStack]==0){
      TopStack--;
      if(TopStack<0)break;
    }
    this->nLooseNodesCoeff=TopStack+1;//move top of stack
  }
  }else{
    #pragma omp atomic
    (this->nLooseNodesCoeff)--;
  }
  */
}

/** SerialTree destructor. */
template<int D>
SerialTree<D>::~SerialTree() {
    for (int i = 0; i < this->genNodeCoeffChunks.size(); i++) delete[] this->genNodeCoeffChunks[i];
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[] (char*)(this->nodeChunks[i]);
    for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];
    for (int i = 0; i < this->genNodeChunks.size(); i++) delete[] (char*)(this->genNodeChunks[i]);

    delete[] this->nodeStackStatus;
    delete[] this->genNodeStackStatus;
    delete[] this->looseNodeStackStatus;

    delete[] this->looseNodeCoeff;
    delete[] this->looseCoeffStack;
}

template class SerialTree<1>;
template class SerialTree<2>;
template class SerialTree<3>;
