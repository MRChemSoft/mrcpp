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
  * Serial tree nodes are not using the destructors, but explicitely call to DeAllocNodes+DeAllocCoeff or DeAllocGenNodes+DeAllocGenCoeff
  * Gen nodes and loose nodes are not counted with MWTree->[in/de]crementNodeCount()
*/
template<int D>
SerialTree<D>::SerialTree(MWTree<D>* Tree,
                                int max_nodes)
        : nNodes(0),
	  nGenNodes(0),
	  nNodesCoeff(0),
	  nGenNodesCoeff(0),
          lastNode(0),
          lastGenNode(0),
          mwTree_p(0),
          maxNodes(max_nodes),
          maxGenNodes(max_nodes),
          maxNodesCoeff(max_nodes),//to be thought of
          maxLooseNodesCoeff(64),//to be thought of
          maxGenNodesCoeff(max_nodes)//to be thought of
          {

    if(max_nodes==0)maxLooseNodesCoeff = 0;

    //Size for GenNodes chunks. ProjectedNodes will be 8 times larger
    int Sizeperchunk = 1024*1024;// 1 MB small for no waisting place, but large enough so that latency and overhead work is negligible

    this->sizeGenNodeCoeff = (MathUtils::ipow(Tree->getOrder()+1,D));//One block
    this->sizeNodeCoeff =(1<<D)*this->sizeGenNodeCoeff;//TDim  blocks
    println(10, "SizeNode Coeff (kB) " << this->sizeNodeCoeff*sizeof(double)/1024);
    println(10, "SizeGenNode Coeff (kB) " << this->sizeGenNodeCoeff*sizeof(double)/1024);


    maxNodesPerChunk = Sizeperchunk/this->sizeGenNodeCoeff;
    println(10, " max_nodes = " <<max_nodes<<", nodes per chunk = "<<maxNodesPerChunk);


    //indicate occupation of nodes
    this->NodeStackStatus = new int[this->maxNodes+1];

    this->LooseNodeCoeff = new double[this->maxLooseNodesCoeff*this->sizeNodeCoeff];

    //  this->SGenData = new double[(this->maxGenNodes*sizeGenNode)/sizeof(double)];
    //indicate occupation of Gen nodes
    this->GenNodeStackStatus = new int[this->maxGenNodes+1];

    //LooseNodeCoeff pointers
    this->LooseCoeffStack = new double * [this->maxLooseNodesCoeff];
    //indicate occupied LooseNodeCoeff
    this->LooseCoeffStackStatus = new int[this->maxLooseNodesCoeff+1];


    //some useful pointers to positions in SData
    //this->lastNode = (ProjectedNode<D>*) (this->SData + this->sizeTreeMeta/sizeof(double));
    this->lastNode = (ProjectedNode<D>*) this->SNodes;
    this->firstNode = (double*)this->lastNode;//constant start of node data
    //this->lastNodeCoeff = (double*) (this->SData+this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNodeMeta/sizeof(double));//start after the metadata
    this->lastNodeCoeff = this->SNodesCoeff;//Start of coefficients
    this->firstNodeCoeff = this->lastNodeCoeff;//constant start of coeff data

    //    this->lastGenNode = (GenNode<D>*) (this->SGenData);
    //this->lastGenNodeCoeff = this->SGenData + (this->maxGenNodes*sizeGenNodeMeta)/sizeof(double);//start at after GenNode Meta part
    this->lastGenNode = this->SGenNodes;
    this->lastGenNodeCoeff = this->SGenNodesCoeff;

    this->nNodesCoeff=-1;//add 1 before each allocation
    this->nLooseNodesCoeff=-1;//add 1 before each allocation
    this->nGenNodesCoeff=-1;//add 1 before each allocation

    //initialize stacks
    for (int i = 0; i <maxNodes;i++){
      this->NodeStackStatus[i] = 0;//0=unoccupied
    }
    this->NodeStackStatus[maxNodes] = -1;//=unavailable


     //initialize stacks
    for (int i = 0; i <this->maxGenNodes;i++){
      this->GenNodeStackStatus[i] = 0;//0=unoccupied
    }
    this->GenNodeStackStatus[this->maxGenNodes] = -1;//=unavailable

    for (int i = 0; i <this->maxLooseNodesCoeff;i++){
      this->LooseCoeffStack[i] = this->LooseNodeCoeff+i*this->sizeNodeCoeff;
      this->LooseCoeffStackStatus[i] = 0;//0=unoccupied
    }
    this->LooseCoeffStackStatus[this->maxLooseNodesCoeff]=-1;//-1=unavailable

    this->mwTree_p = Tree;//pointer to parent tree
    this->mwTree_p->serialTree_p = this;//must be defined before nodes are created

#ifdef HAVE_OPENMP
    omp_init_lock(&Stree_lock);
#endif

    //make virtual table pointers
    FunctionTree<D>* FTree = static_cast<FunctionTree<D>*> (this->mwTree_p);
    const NodeIndex<D> &nIdxtmp = Tree->getRootBox().getNodeIndex(0);
    ProjectedNode<D>* tempnode = new ProjectedNode<D>(*FTree, nIdxtmp);
    this->cvptr_ProjectedNode =  *(char**)(tempnode);
    GenNode<D>* tempGenNode = new GenNode<D>(*tempnode, 0);
    this->cvptr_ProjectedNode =  *(char**)(tempnode);
    this->cvptr_GenNode =  *(char**)(tempGenNode);

    delete tempnode; 
    delete tempGenNode; 
    //alloc root nodes
    MWNode<D> **roots = Tree->getRootBox().getNodes();
    int NodeIx;
    for (int rIdx = 0; rIdx < Tree->getRootBox().size(); rIdx++) {
        const NodeIndex<D> &nIdx = Tree->getRootBox().getNodeIndex(rIdx);
        roots[rIdx] = createSnode(nIdx);
	roots[rIdx]->setIsRootNode(); 
	roots[rIdx]->clearIsGenNode(); 
	println(10, rIdx<<" root idx "<<roots[rIdx]->nodeIndex);
	//	roots[rIdx] = new (allocNodes(1, &NodeIx)) ProjectedNode<D>(*this->getTree(), nIdx);
	//roots[rIdx]->SNodeIx = NodeIx;
    }
    println(10, "tree at " << Tree<<" first root at "<<(Tree->getRootBox().getNodes())[0]<<" SerialTree at "<<this);
    //cout<<MPI_rank<<" root 0 cvptr " << (double*)(*(char**)(roots[0]))<<endl;

    Tree->resetEndNodeTable();


}
/** Overwrite all pointers defined in the tree.
  * Necessary after sending the tree 
  * could be optimized. Should reset other counters? (GenNodes...) */
template<int D>
void SerialTree<D>::RewritePointers(int* &STreeMeta){
  int DepthMax = 100;
  MWNode<D>* stack[DepthMax*8];
  int  slen = 0, counter = 0;

  this->nNodes = STreeMeta[0];
  this->nNodesCoeff = STreeMeta[1];
  this->nGenNodes = 0;
  this->nGenNodesCoeff = -1;
  this->nLooseNodesCoeff = -1;

  //  ptrdiff_t d_p_shift = this->SData - *((double**)this->SData);//actual - written adress
  ptrdiff_t d_p_shift = 0;//actual - written adress of double
  ptrdiff_t n_p_shift = 0;//actual - written adress of node
  ptrdiff_t t_p_shift = 0;//actual - written adress of tree
  cout<<"start of data at "<<this->firstNode<<" pointer at "<<" diff "<<d_p_shift<<endl;
  cout<<" size tree "<< this->getTree()->nNodes <<" tag "<<this<<endl;

  //reinitialize stacks
  for (int i = 0; i <this->maxNodes;i++){
    this->NodeStackStatus[i] = 0;
  }
  for (int i = 0; i <this->maxNodesCoeff;i++){
    this->CoeffStackStatus[i] = 0;//0=unoccupied
  }
  this->CoeffStackStatus[this->maxNodesCoeff]=-1;//-1=unavailable

  for (int i = 0; i <this->maxGenNodes;i++){
    this->GenNodeStackStatus[i] = 0;//0=unoccupied
  }
  this->GenNodeStackStatus[this->maxGenNodes] = -1;//=unavailable

  for (int i = 0; i <this->maxGenNodesCoeff;i++){
    this->GenCoeffStackStatus[i] = 0;//0=unoccupied
  }
  this->GenCoeffStackStatus[this->maxGenNodesCoeff] = -1;//-1=unavailable

  for (int i = 0; i <this->maxLooseNodesCoeff;i++){
    this->LooseCoeffStackStatus[i] = 0;//0=unoccupied
  }
  this->LooseCoeffStackStatus[this->maxLooseNodesCoeff]=-1;//-1=unavailable

  //  for (int i = 0; i <= this->nNodesCoeff; i++){
  //  CoeffStack[i] += d_p_shift;
  //}

  FunctionTree<D>* Tree = static_cast<FunctionTree<D>*> (this->mwTree_p);

  NodeBox<D> &rBox = Tree->getRootBox();
  MWNode<D> **roots = rBox.getNodes();
  //const NodeIndex<D> &nIdx = Tree->getRootBox().getNodeIndex(0);

  d_p_shift = (this->firstNodeCoeff)-(roots[0]->coefs);
  n_p_shift =  ((MWNode<D> *) ( this->firstNodeCoeff)) - ((MWNode<D> *)(roots[0]->coefs));
  //n_p_shift =  (d_p_shift*sizeof(double))/sizeof(MWNode<D>);//be careful with ptrdiff_t arithmetic!
  t_p_shift =  ((MWTree<D> *) ( this->firstNodeCoeff)) - ((MWTree<D> *)(roots[0]->coefs));
  cout<<"d_p_shift "<<d_p_shift<<" as tree* "<<t_p_shift<<" as node* "<<n_p_shift<<" tag "<<this<<endl;
  //test consistency: could use explicit shift instead
  //if((n_p_shift*sizeof(MWNode<D>))/sizeof(double) != d_p_shift)MSG_FATAL("Serial Tree: wrong n_p_shift");
  //if((t_p_shift*sizeof(MWTree<D>))/sizeof(double) != d_p_shift)cout<<" wrong t_p_shift "<<endl;


  for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
   stack[slen++] =  roots[rIdx];//roots address are in MWtree which is not overwritten
  }
  this->getTree()->nNodes = 0;
  while (slen) {
      this->getTree()->nNodes++;
      MWNode<D>* fpos = stack[--slen];
      if(fpos->getNChildren()){
 	for (int i = 0; i < fpos->getNChildren(); i++) {
	  fpos->children[i] = (MWNode<D>*)(((double*)fpos->children[i]) + d_p_shift);
	  stack[slen++] = fpos->children[i];
	}
      }
      //set virtual table. Assumes a Projected Node!
      //cout<<"cvptr BEFORE "<<(double*)(*(char**)(fpos))<<" AFTER "<<(double*)cvptr<<endl;
      *(char**)(fpos) = this->cvptr_ProjectedNode;

      fpos->parent = (MWNode<D>*)(((double*)fpos->parent) + d_p_shift);// += n_p_shift not safe!
      //fpos->tree += t_p_shift;//NB: does not work!//(MWTree<D>*)(((double*)fpos->tree) + d_p_shift);
      fpos->tree = this->getTree();
      fpos->coefs = this->firstNodeCoeff+(fpos->SNodeIx)*this->sizeNodeCoeff;
      this->NodeStackStatus[fpos->SNodeIx] = 1;

   }
  this->getTree()->resetEndNodeTable();

}
/** Adds two trees.
  * Generate missing nodes on the flight and adds all nodes */
/*template<int D>
void SerialTree<D>::SerialTreeAdd(double c, FunctionTree<D>* &TreeB, FunctionTree<D>* &TreeC){
  println(0, " SerialTreeAdd ");
    int DepthMax = 100;
    int  slen = 0, counter = 0;
    int  slenA = 0, counterA = 0, slenB = 0, counterB = 0;
    MWNode<D>* stackA[DepthMax*8];
    MWNode<D>* stackB[DepthMax*8];
    FunctionTree<D>* TreeA=this->getTree();
    int N_GenCoeff=TreeA->getKp1_d();
    int N_Coeff=N_GenCoeff*TreeA->getTDim();
    int Children_Stride = this->sizeNodeCoeff;
    double Tsum=0.0;

    if(TreeA->getRootBox().size()!= TreeB->getRootBox().size())MSG_FATAL("Number of root nodes must be equal for now");
 
    //put root nodes on stack
    for (int rIdxA = 0; rIdxA < TreeA->getRootBox().size(); rIdxA++) stackA[slenA++] = TreeA->getRootBox().getNodes()[rIdxA];
    for (int rIdxB = 0; rIdxB < TreeB->getRootBox().size(); rIdxB++) stackB[slenB++] = TreeB->getRootBox().getNodes()[rIdxB];
 
    Timer timer,t0,t1,t2,t3;
    double* cA;
    double* cB;
    double* cC;
    timer.start();
    t1.start();
    t2.start();
    Tsum=0.0;
    counter=0;
    counterA=0;
    counterB=0;
    while (slenA) {
      counter++;
      MWNode<D>* fposA = stackA[--slenA];
      MWNode<D>* fposB = stackB[slenA];
      //	  println(0, " treating   " <<stackA[slenA]->getRank()<<" with slenA "<<slenA);
      if(fposA->getNChildren()+fposB->getNChildren()){
	cA=&(fposA->getCoefs()(0));
	if(fposA->getNChildren()==0){
	  t1.resume();
	  GenS_nodes(fposA);
	  t1.stop();
	  counterB++;
	}
	if(fposB->getNChildren()==0){
	  t1.resume();
	  cB=&(fposB->getCoefs()(0));
	  GenS_nodes(fposB);
	  counterB++;
	  t1.stop();
	}
	for (int i = 0; i < fposA->getNChildren(); i++) {
	  stackA[slenA] = fposA->children[i];
	  stackB[slenA++] = fposB->children[i];
	}
      }

      if(1){
	counterA++;
	cA=&(fposA->getCoefs()(0));
	cB=&(fposB->getCoefs()(0));
	//virtual table test
	//char* cvptr =  *(char**)(fposA);
	//println(0, fposA->getRank()<<" virtual pointer "<<(double*) cvptr <<" "<<fposA->NodeCoeffIx<<" "<<fposA->GenNodeCoeffIx);  
	
	t2.resume();
	if(fposA->hasWCoefs()){
	  if(fposB->hasWCoefs()){for (int i=0; i<N_Coeff; i++)cA[i]+=c*cB[i];
	  }else{for (int i=0; i<N_GenCoeff; i++)cA[i]+=c*cB[i];
	  }	  	
	}else{
	  for (int i=0; i<N_GenCoeff; i++)cA[i]+=c*cB[i];
	  if(fposB->hasWCoefs()){for (int i=N_GenCoeff; i<N_Coeff; i++)cA[i]=c*cB[i];
	  }else{println(0, slenA<<" Adding two generated nodes? " );//should not happen	  
	  }	  
	}
	fposA->setHasWCoefs();
	fposA->calcNorms();
	if(fposA->getNChildren()==0)Tsum+=fposA->getSquareNorm();
	t2.stop();
      }
    }
    println(0, " summed "<<counterA<<" generated "<<counterB<<" looped "<<counter);
    println(0, " squarenorm "<<Tsum);


    this->getTree()->squareNorm=Tsum;
    println(0, " time generate     " << t1);
    println(0, " time add coef     " << t2);
    timer.stop();
    println(0, " time Sadd     " << timer);

    this->getTree()->resetEndNodeTable();
    cout<<"sending TreeAB with Nnodes "<<this->nNodes<<endl;
#ifdef HAVE_MPI
    //    if(MPI_size == 2)SendRcv_SerialTree(TreeA, 0, 1, 44, MPI_COMM_WORLD);
#endif


}*/

/** Adds two trees.
  * Generate missing nodes on the flight and "compress" the ancestor  from summed nodes also on the flight */
/*template<int D>
void SerialTree<D>::SerialTreeAdd_Up(double c, FunctionTree<D>* &TreeB, FunctionTree<D>* &TreeC){
  println(0, " SerialTreeAdd ");
  //to do: traverse from low to high rank
    int DepthMax = 100;
    int  slen = 0, counter = 0;
    int  slenA = 0, counterA = 0;
    MWNode<D>* stackA[DepthMax*8];
    int  slenB = 0, counterB = 0;
    MWNode<D>* stackB[DepthMax*8];
    FunctionTree<D>* TreeA=this->getTree();
    int N_GenCoeff=TreeA->getKp1_d();
    int tDim=TreeA->getTDim();
    int N_Coeff=N_GenCoeff*tDim;
    //put root nodes on stack
    NodeBox<D> &rBoxA = TreeA->getRootBox();
    MWNode<D> **rootsA = rBoxA.getNodes();
    NodeBox<D> &rBoxB = TreeB->getRootBox();
    MWNode<D> **rootsB = rBoxB.getNodes();
    int Children_Stride = this->sizeNodeCoeff;
    MWNode<D>* fposAA;
    MWNode<D>* fposBB;

    int flag=1;
    double Tsum=0.0;

   slenA=-1;
   counterA=0;
   slenB=-1;
   counterB=0;
    if(rBoxB.size()!=rBoxA.size())MSG_FATAL("Number of root nodes must be equal for now");
    for (int rIdxA = 0; rIdxA < rBoxA.size(); rIdxA++) {
       const NodeIndex<D> &nIdx = rBoxA.getNodeIndex(rIdxA);
	stackA[++slenA] = TreeA->findNode(nIdx);
    }
    for (int rIdxB = 0; rIdxB < rBoxB.size(); rIdxB++) {
        const NodeIndex<D> &nIdx = rBoxB.getNodeIndex(rIdxB);
	stackB[++slenB] = TreeB->findNode(nIdx);
    }

    Timer timer,t0,t1,t2,t3;
    double* cA;
    double* cB;
    double* cC;
    bool downwards = true;//travers in direction children
    timer.resume();
    t2.resume();
    Tsum=0.0;
    counter=0;
    counterA=0;
    counterB=0;
    while (slenA>=0) {
      counter++;
      MWNode<D>* fposA = stackA[slenA];
      MWNode<D>* fposB = stackB[slenA];
      //	println(0, slenA<<" TreeA children   " <<fposA->getNChildren()<<" TreeB children   " <<fposB->getNChildren());

      //println(0, " Treating   " <<fposA->getRank()<<" counter "<<counter);
      
      if(fposA->getNChildren()+fposB->getNChildren() and downwards){
	//println(0, fposA->getRank()<<" Tree children   " <<fposA->getNChildren()<<" TreeB children   " <<fposB->getNChildren()<<" counter "<<counter);
	cA=&(fposA->getCoefs()(0));
	//	  println(0, "NodeA orig at  "<< fposA<<" coeff at "<<cA);
	if(fposA->getNChildren()==0){
	  t1.resume();
	  GenS_nodes(fposA);
	  t1.stop();
	}
	if(fposB->getNChildren()==0){
	  //println(0, slenB<<" BgenChildren()  ");
	  t1.resume();
	  cB=&(fposB->getCoefs()(0));
	  GenS_nodes(fposB);

	  if(flag){
	    flag=0;
	    cB=&(fposB->children[0]->getCoefs()(0));
	  }
	  t1.stop();

	  //println(0, slenB<<" BgenChildren() succesful ");
	}
	//	println(0, slenA<<" TreeA newchildren   " <<fposA->getNChildren()<<" TreeB newchildren   " <<fposB->getNChildren());
	for (int i = 0; i < fposA->getNChildren(); i++) {
	  stackA[++slenA] = fposA->children[i];
	  stackB[slenA] = fposB->children[i];
	  //	  println(0, " child   " <<stackA[slenA]->getRank()<<" with slenA "<<slenA);
	  downwards = true;
	}
	//println(0, " oldest   " <<stackA[slenA]->getRank()<<" with slenA "<<slenA);
      }else{
	bool youngestchild = false;
	if(fposA->parent!=0) youngestchild = (fposA->parent->children[0]==fposA);
	//      println(0, slenA<<"  youngestchild = "<<youngestchild);
	if(youngestchild or slenA==0){
	  //println(0, slenA<<" adding siblings of  "<<fposA->getRank());
	  //all siblings are ready
	  int siblingsize=tDim;
	  if(slenA==0)siblingsize=rBoxA.size();
	  for (int ichild = 0; ichild < siblingsize; ichild++) {
	    if(fposA->parent!=0){
	      fposAA=fposA->parent->children[ichild];
	      fposBB=fposB->parent->children[ichild];
	    }else{
	      fposAA=TreeA->findNode(rBoxA.getNodeIndex(ichild));
	      fposBB=TreeB->findNode(rBoxB.getNodeIndex(ichild));
	    }
	    if(fposAA->getNChildren()==0){
	      //sum all siblings
	      //      if(fposA->getNChildren() + fposB->getNChildren() ==0){
	      //if(1){
	      cA=&(fposAA->getCoefs()(0));
	      cB=&(fposBB->getCoefs()(0));
	      	      
	      t2.resume();
	      if(fposAA->isGenNode()){
		
		//for (int i=0; i<N_GenCoeff; i++)fposA->getCoefs()(i)+=c*fposB->getCoefs()(i);
		for (int i=0; i<N_GenCoeff; i++)cA[i]+=c*cB[i];
		if(fposBB->isGenNode()){
		  //should not happen	  
		  println(0, slenA<<" Adding two generated nodes? " );
		  //for (int i=N_GenCoeff; i<N_Coeff; i++)fposA->coefs[i]=0.0;
		}else{
		  for (int i=N_GenCoeff; i<N_Coeff; i++)cA[i]=c*cB[i];
		}	  
	      }else{
		if(fposBB->isGenNode()){
		  for (int i=0; i<N_GenCoeff; i++)cA[i]+=c*cB[i];
		}else{
		  for (int i=0; i<N_Coeff; i++)cA[i]+=c*cB[i];
		}	  	
	      }
	      t2.stop();
	      fposAA->calcNorms();
	      Tsum+=fposAA->getSquareNorm();
	      println(0, " rank   " <<fposAA->getRank()<<" norm  "<<fposAA->getSquareNorm());
	      counterA++;
	    }
	  }
	  if(slenA>0){
	    //make parent
	    t3.resume();
	    S_mwTransformBack(&(fposA->getCoefs()(0)), &(fposA->parent->getCoefs()(0)), Children_Stride); 
	    fposA->parent->calcNorms();
	    t3.stop();
	    //println(0, " made parent "<<fposA->parent->getRank());
	    counterB++;
	  }
	  downwards = false;
	  slenA--;
	}else{
	  //println(0, " did nothing for "<<stackA[slenA]->getRank());
	  slenA--;
	  downwards = true;
	}
      }
    }
    println(0, " summed "<<counterA<<" generated "<<counterB<<" looped "<<counter);
    println(0, " squarenorm "<<Tsum);


    this->getTree()->squareNorm=Tsum;
    println(0, " time generate     " << t1);
    println(0, " time add coef     " << t2);
    t3.resume();
    //this->S_mwTreeTransformUp();
     //this->getTree()->mwTransform(BottomUp);
    t3.stop();
    println(0, " time TransformUp    " << t3);
    timer.stop();
    println(0, " time Sadd     " << timer);

   println(0, "TreeAB Nodes   " <<this->nNodes<<" squarenorm "<<Tsum);

   }*/

/** Make 8 children nodes with scaling coefficients from parent
 * Does not put 0 on wavelets
 */
template<int D>
void SerialTree<D>::GenS_nodes(MWNode<D>* Node){

  bool ReadOnlyScalingCoeff=true;

  //  if(Node->isGenNode())ReadOnlyScalingCoeff=true;
  if(Node->hasWCoefs())ReadOnlyScalingCoeff=false;
  double* cA;

  Node->genChildren();//will make children and allocate coeffs, but without setting values for coeffs.


  double* coeffin  = Node->coefs;
  double* coeffout = Node->children[0]->coefs;

  int Children_Stride = this->sizeGenNodeCoeff;
  S_mwTransform(coeffin, coeffout, ReadOnlyScalingCoeff, Children_Stride);

}


/** Make children scaling coefficients from parent
 * Other node info are not used/set
 * this routine works only for D=3.
 * coeff_in is not modified.
 * The output is written directly into the 8 children scaling coefficients. 
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
template<int D>
void SerialTree<D>::S_mwTransform(double* coeff_in, double* coeff_out, bool ReadOnlyScalingCoeff, int Children_Stride, bool b_overwrite) {

  int operation = Reconstruction;
  int kp1 = this->getTree()->getKp1();
  int tDim = (1<<D);
  int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
  int kp1_d = kp1_dm1*kp1;//
  const MWFilter &filter = this->getTree()->getMRA().getFilter();
  //VectorXd &result = Parent->tree->getTmpMWCoefs();
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
      /* Operate in direction i only if the bits along other
       * directions are identical. The bit of the direction we
       * operate on determines the appropriate filter/operator */
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
      /* Operate in direction i only if the bits along other
       * directions are identical. The bit of the direction we
       * operate on determines the appropriate filter/operator */
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
	/* Operate in direction i only if the bits along other
	 * directions are identical. The bit of the direction we
	 * operate on determines the appropriate filter/operator */
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
}

/** Regenerate all s/d-coeffs by backtransformation, starting at the bottom and
  * thus purifying all coefficients.
  * not yet fully optimized.
 */
/*template<int D>
void SerialTree<D>::S_mwTreeTransformUp() {
    Timer t0,t1,t2,t3;
    int DepthMax = 100;
    MWTree<D>* Tree=this->getTree();
    int status[this->nNodes];
    MWNode<D>* stack[DepthMax*8];
    int  slen = 0, counter = 0, counter2 = 0;
    NodeBox<D> &rBox = Tree->getRootBox();
    MWNode<D> **roots = rBox.getNodes();
    int Children_Stride = this->sizeGenNodeCoeff;
    for (int i = 0; i < this->nNodes; i++) status[i]=0;
    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
       const NodeIndex<D> &nIdx = rBox.getNodeIndex(rIdx);
	stack[++slen] = Tree->findNode(nIdx);
	MWNode<D>* fpos = stack[slen];
	if(fpos->getNChildren()>0){
	  status[fpos->getRank()]=0;
	}else {
	  status[fpos->getRank()]=1;//finished
	}
    }
    while (slen) {
      counter++;
	MWNode<D>* fpos = stack[slen];
	if(fpos->getNChildren()>0 and status[fpos->getRank()]==0){
	  int childok=0;
	  for (int i = 0; i < fpos->getNChildren(); i++) {
	    if(status[fpos->children[i]->getRank()]>0 or fpos->children[i]->getNChildren()==0){
	      childok++;
	    }else{
	      if(fpos->children[i]->getNChildren()>0){
		stack[++slen] = fpos->children[i];
	      }
	    }
	  }
	  if(childok==fpos->getNChildren()){
	    //Ready to compress fpos from children
	    t0.resume();
	    S_mwTransformBack(&(fpos->children[0]->getCoefs()(0)), &(fpos->getCoefs()(0)), Children_Stride); 
	    t0.stop();
	    status[fpos->getRank()]=1;//finished     
	    counter2++;
	  }
	    
	}else {
	  slen--;
	  status[fpos->getRank()]=1;	  
	}
    }
    println(0, " time   S_mwTransformBack   " << t0);
    println(0, counter2<<" nodes recompressed, out of "<<this->nNodes);

    
    }*/

/** Make parent from children scaling coefficients
 * Other node info are not used/set
 * coeff_in is not modified.
 * The output is read directly from the 8 children scaling coefficients. 
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
/*template<int D>
void SerialTree<D>::S_mwTransformBack(double* coeff_in, double* coeff_out, int Children_Stride) {
  if(D!=3)MSG_FATAL("S_mwtransform Only D=3 implemented now!");
  int operation = Compression;
  int kp1 = this->getTree()->getKp1();
  int tDim = (1<<D);
  int kp1_dm1 = kp1*kp1;//MathUtils::ipow(Parent->getKp1(), D - 1)
  int kp1_d = kp1_dm1*kp1;//
  const MWFilter &filter = this->getTree()->getMRA().getFilter();
  //VectorXd &result = Parent->tree->getTmpMWCoefs();
  double overwrite = 0.0;
  double *tmp;

  double tmpcoeff[kp1_d*tDim];
  //double tmp2coeff[kp1_d*tDim];
  

  int ftlim=tDim;
  int ftlim2=tDim;
  int ftlim3=tDim;

  int i = 0;
  int mask = 1;
  for (int gt = 0; gt < tDim; gt++) {
    double *out = coeff_out + gt * kp1_d;
    for (int ft = 0; ft < ftlim; ft++) {
      /* Operate in direction i only if the bits along other
       * directions are identical. The bit of the direction we
       * operate on determines the appropriate filter/operator */
/*      if ((gt | mask) == (ft | mask)) {
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
      /* Operate in direction i only if the bits along other
       * directions are identical. The bit of the direction we
       * operate on determines the appropriate filter/operator */
/*     if ((gt | mask) == (ft | mask)) {
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
      /* Operate in direction i only if the bits along other
       * directions are identical. The bit of the direction we
       * operate on determines the appropriate filter/operator */
/*     if ((gt | mask) == (ft | mask)) {
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

  }*/


/** Allocating a Projected Serial (root) node.
  *
  * This routine creates an empty ProjectedNode node
  * with the appropriate translation */
template<int D>
ProjectedNode<D>* SerialTree<D>::createSnode(const NodeIndex<D> & nIdx) {

  int NodeIx;
  double *coefs_p;
  ProjectedNode<D>* newNode=this->allocNodes(1, &NodeIx, &coefs_p);

  *(char**)(newNode) = this->cvptr_ProjectedNode;
  newNode->SNodeIx = NodeIx;
  newNode->tree = this->mwTree_p;
  newNode->parent = 0;
  newNode->nodeIndex = nIdx;
  newNode->hilbertPath = HilbertPath<D>();
  newNode->squareNorm = -1.0;
  newNode->status = 0;
  newNode->n_coefs = 0;
  newNode->coefs = 0;
  newNode->clearNorms();
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

}


//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* SerialTree<D>::allocNodes(int nAlloc, int* NodeIx, double** coefs_p) {

  *NodeIx = this->nNodes;
  int ChunkIx = *NodeIx%(this->maxNodesPerChunk);

  if(ChunkIx == 0  or  ChunkIx+nAlloc > maxNodesPerChunk ){

    if (this->nNodes+nAlloc >= this->maxNodes){
      println(0, "maxNodes exceeded " << this->maxNodes);
      MSG_FATAL("maxNodes exceeded ");
    } 

    //we want nodes allocated simultaneously to be allocated on the same pice.
    //possibly jump over the last nodes from the old chunk
    this->nNodes=this->maxNodesPerChunk*((this->nNodes+nAlloc-1)/maxNodesPerChunk);//start of next chunk

    int chunk = this->nNodes/maxNodesPerChunk;//find the right chunk

    //careful: NodeChunks.size() is an unsigned int
    if(chunk+1 > this->NodeChunks.size()){
	//need to allocate new chunk
	this->SNodes = (ProjectedNode<D>*) new char[this->maxNodesPerChunk*sizeof(ProjectedNode<D>)];
	this->NodeChunks.push_back(this->SNodes);
        this->SNodesCoeff = new double[this->sizeNodeCoeff*this->maxNodesPerChunk];
        NodeCoeffChunks.push_back(SNodesCoeff);
	if(chunk%100==99 and D==3)println(10,endl<<" number of nodes "<<this->nNodes <<",number of Nodechunks now " << this->NodeChunks.size()<<", total size coeff  (MB) "<<(this->nNodes/1024) * this->sizeNodeCoeff/128);
      }
      this->lastNode = this->NodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);	
      *NodeIx = this->nNodes;
      ChunkIx = *NodeIx%(this->maxNodesPerChunk);

    }
    assert((this->nNodes+nAlloc-1)/maxNodesPerChunk < this->NodeChunks.size());

    ProjectedNode<D>* newNode  = this->lastNode;
    ProjectedNode<D>* newNode_cp  = newNode;
    *coefs_p = this->SNodesCoeff + ChunkIx*this->sizeNodeCoeff;
 
    //We can have sizeNodeMeta with different size than ProjectedNode
    //println(0, "new size meta " << this->nNodes);
    for (int i = 0; i<nAlloc; i++){
      if (this->NodeStackStatus[*NodeIx+i]!=0)
	println(0, *NodeIx+i<<" NodeStackStatus: not available " << this->NodeStackStatus[*NodeIx+i]);
      this->NodeStackStatus[*NodeIx+i]=1;
      newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->lastNode+=nAlloc;
    return newNode;
}

template<int D>
void SerialTree<D>::DeAllocNodes(int SNodeIx) {
  if (this->nNodes <0){
    println(0, "minNodes exceeded " << this->nNodes);
    this->nNodes++;
  }
  this->NodeStackStatus[SNodeIx]=0;//mark as available
  if(SNodeIx==this->nNodes-1){//top of stack
    int TopStack=this->nNodes;
    while(this->NodeStackStatus[TopStack-1]==0){
      TopStack--;
      if(TopStack<1)break;
    }
    this->nNodes=TopStack;//move top of stack
    //has to redefine lastGenNode                                                                                                                   
    int chunk = this->nNodes/maxNodesPerChunk;//find the right chunk                                                                               
    this->lastNode = this->NodeChunks[chunk] + this->nNodes%(this->maxNodesPerChunk);

  }
 }

//return pointer to the last active node or NULL if failed
template<int D>
GenNode<D>* SerialTree<D>::allocGenNodes(int nAlloc, int* NodeIx, double** coefs_p) {

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
}

template<int D>
void SerialTree<D>::DeAllocGenNodes(int SNodeIx) {
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
 }

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocCoeff(int nAllocCoeff, MWNode<D>* node) {

  if(node->isLooseNode()){
    return this->allocLooseCoeff(nAllocCoeff, node);
  }else{
     println(0, "ERROR SHOUDL NOT ALLOCATE COEFF HERE " );
    omp_set_lock(&Stree_lock);
    if (nAllocCoeff!=1<<D) MSG_FATAL("Only 2**D implemented now!");
    if (this->nNodesCoeff+1 >= this->maxNodesCoeff ){
      println(0, "maxNodesCoeff exceeded " << this->maxNodesCoeff);
      MSG_FATAL("maxNodesCoeff exceeded ");
      omp_unset_lock(&Stree_lock);
      return 0;
    } else if( CoeffStackStatus[this->nNodesCoeff+1]!=0){
      println(0, this->nNodesCoeff+1<<" CoeffStackStatus: not available " <<CoeffStackStatus[this->nNodesCoeff+1] <<" tree at "<<this->mwTree_p);      
      omp_unset_lock(&Stree_lock);
      return 0;
    }else{

    
    if((this->nNodesCoeff)/maxNodesPerChunk != (this->nNodesCoeff+1)/maxNodesPerChunk){
      //we want nodes allocated simultaneously to be allocated on the same chunk.
       //jump over the last nodes from the old chunk
      this->nNodesCoeff=this->maxNodesPerChunk*((this->nNodesCoeff+1)/maxNodesPerChunk)-1;//start of next chunk
    }
    if((this->nNodesCoeff+1)/maxNodesPerChunk >= NodeCoeffChunks.size()){
      //need to allocate new chunk
      this->SNodesCoeff = new double[this->sizeNodeCoeff*this->maxNodesPerChunk];

      NodeCoeffChunks.push_back(SNodesCoeff);
      println(0, "allocated new chunk for Snodes Coeff. number of chunks now " << this->NodeCoeffChunks.size());
    }
    this->nNodesCoeff += 1;//nAllocCoeff;

      double* NodeCoeff=0;//this->CoeffStack[this->nNodesCoeff];
      this->CoeffStackStatus[this->nNodesCoeff]=1;
      node->SNodeIx = this->nNodesCoeff;
      //node->NodeCoeffIx = node->SNodeIx;
      omp_unset_lock(&Stree_lock);
      return NodeCoeff;
    }
  }
}

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocLooseCoeff(int nAllocCoeff, MWNode<D>* node) {

  if (nAllocCoeff!=1<<D) MSG_FATAL("Only 2**D implemented now!");

  //Each omp thread use own part of array, no locks!

#pragma omp atomic
  (this->nLooseNodesCoeff)++;

  if(this->nLooseNodesCoeff>8){
    cout<<" Loose nodes now: "<<this->nLooseNodesCoeff<<endl;
  }
  if (this->nLooseNodesCoeff >= this->maxLooseNodesCoeff ){
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
    return LooseNodeCoeff;
  }

}


//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocCoeff(int Index) {

    MSG_FATAL("NodesCoeff should be allocated with nodes ");
  omp_set_lock(&Stree_lock);
  this->nNodesCoeff += 1;
  if (this->nNodesCoeff >= this->maxNodesCoeff ){
    println(0, "maxNodesCoeff exceeded " << this->maxNodesCoeff);
    MSG_FATAL("maxNodesCoeff exceeded ");
    omp_unset_lock(&Stree_lock);
    return 0;
  } else if( CoeffStackStatus[Index]!=0){
    println(0, Index<<" CoeffStackStatus: not available " <<CoeffStackStatus[Index] );      
    MSG_FATAL("CoeffStackStatus: not available  ");
    omp_unset_lock(&Stree_lock);
    return 0;
  }else{

 
    if(Index >= this->maxNodesPerChunk*NodeCoeffChunks.size()){
      //need to allocate new chunk
      this->SNodesCoeff = new double[this->sizeNodeCoeff*this->maxNodesPerChunk];
      NodeCoeffChunks.push_back(SNodesCoeff);
      println(0, "allocated new chunk for Snodes Coeff I. number of chunks now " << this->NodeCoeffChunks.size());
    }

    double* NodeCoeff=0;//this->CoeffStack[Index];
    this->CoeffStackStatus[Index]=1;
    omp_unset_lock(&Stree_lock);
    return NodeCoeff;
  }
}

//"Deallocate" the stack
template<int D>
void SerialTree<D>::DeAllocCoeff(int DeallocIx) {  
  omp_set_lock(&Stree_lock);
    println(0, "ERROR SHOULD NOT BE USED " << DeallocIx);
  if (this->CoeffStackStatus[DeallocIx]==0){
    println(0, "deleting already unallocated coeff " << DeallocIx<<" "<<this->LooseCoeffStackStatus[DeallocIx]<<" nNodes "<<this->nNodes<<" nNodesCoeff "<<this->nNodesCoeff<<" nGenNodesCoeff "<<this->nGenNodesCoeff);
  }
  this->CoeffStackStatus[DeallocIx]=0;//mark as available
  if(DeallocIx==this->nNodesCoeff){//top of stack
    int TopStack=this->nNodesCoeff;
    while(this->CoeffStackStatus[TopStack]==0){
      TopStack--;
      if(TopStack<1)break;
    }
    this->nNodesCoeff=TopStack;//move top of stack
  }
  omp_unset_lock(&Stree_lock);
}


//"Deallocate" the stack
template<int D>
void SerialTree<D>::DeAllocLooseCoeff(int DeallocIx) {  
  if (this->LooseCoeffStackStatus[DeallocIx]==0){
    println(0, "deleting already unallocated loose coeff " << DeallocIx);
  }
  this->LooseCoeffStackStatus[DeallocIx]=0;//mark as available
  if(omp_get_num_threads()==1){
  if(DeallocIx==this->nLooseNodesCoeff){//top of stack
    int TopStack=this->nLooseNodesCoeff;
    while(this->LooseCoeffStackStatus[TopStack]==0){
      TopStack--;
      if(TopStack<1)break;
    }
    this->nLooseNodesCoeff=TopStack;//move top of stack
  }
  }else{
    #pragma omp atomic
    (this->nLooseNodesCoeff)--;
  }
}

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocGenCoeff(int nAllocCoeff, MWNode<D>* node) {
  println(0, "ERROR allocating GenNodesCoeff " << node->SNodeIx<<" "<<this);
  if(node->isLooseNode()){
    return this->allocLooseCoeff(nAllocCoeff, node);
  }else{
   omp_set_lock(&Stree_lock);
  if (nAllocCoeff!=1<<D) MSG_FATAL("Only 2**D implemented now!");

  if (this->nGenNodesCoeff+1 >= this->maxGenNodesCoeff ){
    println(0, "maxGenNodesCoeff exceeded " << this->maxGenNodesCoeff);
    MSG_FATAL("maxGenNodesCoeff exceeded ");
    omp_unset_lock(&Stree_lock);
    return 0;
    //  } else if( GenCoeffStackStatus[this->nGenNodesCoeff]!=0){
  } else if( GenCoeffStackStatus[node->SNodeIx]!=0){
    println(0, node->SNodeIx<<" "<<node->SNodeIx<<" A GenCoeffStackStatus: not available " <<GenCoeffStackStatus[node->SNodeIx] <<" tree at "<<this->mwTree_p );      
    omp_unset_lock(&Stree_lock);
    return 0;
  }else{


    if(node->SNodeIx >= this->maxNodesPerChunk*GenNodeCoeffChunks.size()){
      //need to allocate new chunk
      this->SGenNodesCoeff = new double[this->sizeGenNodeCoeff*this->maxNodesPerChunk];
      GenNodeCoeffChunks.push_back(SGenNodesCoeff);
      println(0, "allocated new chunk for SGennodes Coeff. number of chunks now " << this->GenNodeCoeffChunks.size());
    }

    this->nGenNodesCoeff += 1;//nAllocCoeff;


    double* GenNodeCoeff=0;//this->GenCoeffStack[node->SNodeIx];
    this->GenCoeffStackStatus[node->SNodeIx]=1;
    node->SNodeIx = node->SNodeIx;
    //for (int i = 0; i < this->sizeGenNodeCoeff; i++)coeffVectorXd[i]=0.0;
    omp_unset_lock(&Stree_lock);
    return GenNodeCoeff;
  }
  }
}

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocGenCoeff(int Index) {

  println(0, "ERROR allocating GenNodesCoeff Index " << Index);
  omp_set_lock(&Stree_lock);
  this->nGenNodesCoeff += 1;//nAllocCoeff;
  if (this->nGenNodesCoeff >= this->maxGenNodesCoeff ){
    println(0, "maxGenNodesCoeff exceeded " << this->maxGenNodesCoeff);
    MSG_FATAL("maxGenNodesCoeff exceeded ");
    this->nGenNodesCoeff -= 1;//nAllocCoeff;
    omp_unset_lock(&Stree_lock);
    return 0;
    //  } else if( GenCoeffStackStatus[this->nGenNodesCoeff]!=0){
  } else if( GenCoeffStackStatus[Index]!=0){
    println(0, Index<<" GenCoeffStackStatus: not available " <<GenCoeffStackStatus[Index]  <<" tree at "<<this->mwTree_p<<" nGenNodesCoeff "<<this->nGenNodesCoeff<<" nGenNodes "<<this->nGenNodes<<" nNodes "<<this->nNodes);      
    omp_unset_lock(&Stree_lock);
    return 0;
   }else{

    if(Index >= this->maxNodesPerChunk*GenNodeCoeffChunks.size()){
      //need to allocate new chunk
      this->SGenNodesCoeff = new double[this->sizeGenNodeCoeff*this->maxNodesPerChunk];
      GenNodeCoeffChunks.push_back(SGenNodesCoeff);
      println(0, "allocated new chunk for SGennodes Coeff I. number of chunks now " << this->GenNodeCoeffChunks.size());
    }

    double* GenNodeCoeff=0;//this->GenCoeffStack[Index];
    this->GenCoeffStackStatus[Index]=1;
    omp_unset_lock(&Stree_lock);
    return GenNodeCoeff;
  }
}

//"Deallocate" the stack
template<int D>
void SerialTree<D>::DeAllocGenCoeff(int DeallocIx) {
    println(0, "ERROR deleting Gencoeff " << DeallocIx<<" "<<this->mwTree_p);
  omp_set_lock(&Stree_lock);

  if(this->GenCoeffStackStatus[DeallocIx]==0){
    println(0, "deleting already unallocated Gencoeff " << DeallocIx<<" "<<this->mwTree_p);
  }else{
    this->GenCoeffStackStatus[DeallocIx]=0;//mark as available
    if(DeallocIx==this->nGenNodesCoeff){//top of stack
      int TopStack=this->nGenNodesCoeff;
      while(this->CoeffStackStatus[TopStack]==0){
	TopStack--;
	if(TopStack<1)break;
      }
      this->nGenNodesCoeff=TopStack;//move top of stack
    }
  }
  omp_unset_lock(&Stree_lock);
}

/** SerialTree destructor. */
template<int D>
SerialTree<D>::~SerialTree() {
  //println(10, "~SerialTree");
  if(this->maxNodes>0){
    MWNode<D> **roots = this->getTree()->getRootBox().getNodes();
    for (int i = 0; i < this->getTree()->getRootBox().size(); i++) {
        ProjectedNode<D> *node = static_cast<ProjectedNode<D> *>(roots[i]);
	//        node->~ProjectedNode();
        node->deleteChildren();
	this->DeAllocNodes(node->SNodeIx);
	//this->DeAllocCoeff(node->SNodeIx);
	this->mwTree_p->decrementNodeCount(node->getScale());
        roots[i] = 0;
    }

    for (int i = 0; i < this->GenNodeCoeffChunks.size(); i++)delete[] this->GenNodeCoeffChunks[i];
    for (int i = 0; i < this->NodeChunks.size(); i++)delete[] (char*)(this->NodeChunks[i]);
    for (int i = 0; i < this->NodeCoeffChunks.size(); i++)delete[] this->NodeCoeffChunks[i];
    for (int i = 0; i < this->GenNodeChunks.size(); i++)delete[] (char*) (this->GenNodeChunks[i]);

    delete[] this->NodeStackStatus;
    delete[] this->GenNodeStackStatus;

//    delete[] this->CoeffStack;
    //    delete[] this->CoeffStackStatus;

    delete[] this->LooseNodeCoeff;
    delete[] this->LooseCoeffStack;
    delete[] this->LooseCoeffStackStatus;

    //    delete[] this->GenCoeffStack;
    //    delete[] this->GenCoeffStackStatus;

  }else{

    //NB: the pointer values may be wrong, cannot use delete
    //this->SData=0;
    this->SNodes=0;;
    this->SNodesCoeff=0;;
    this->NodeStackStatus=0;

    this->SGenNodes=0;;
    this->SGenNodesCoeff=0;;
    this->GenNodeStackStatus=0;

    this->LooseNodeCoeff=0;
    this->LooseCoeffStack=0;
    this->LooseCoeffStackStatus=0;

  }
}




template class SerialTree<1>;
template class SerialTree<2>;
template class SerialTree<3>;
