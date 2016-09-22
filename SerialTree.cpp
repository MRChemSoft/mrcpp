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

/** SerialTree constructor.
  * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
  * Initializes rootNodes to represent the zero function. */
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
          SData(0),
          maxNodes(max_nodes),
          maxGenNodes(max_nodes),
          maxNodesCoeff(max_nodes),//to be thought of
          maxGenNodesCoeff(max_nodes),//to be thought of
          sizeTreeMeta(0),
          sizeNodeMeta(0),
          sizeGenNodeMeta(0),
          sizeNode(0),
          sizeGenNode(0) {

    println(10, " max_nodes  " <<max_nodes);

    this->sizeNodeMeta =16*((sizeof(ProjectedNode<D>)+15)/16);//we want only multiples of 16
    this->sizeGenNodeMeta =16*((sizeof(GenNode<D>)+15)/16);//we want only multiples of 16
    this->sizeNodeCoeff = 8*(((1<<D)*(MathUtils::ipow(Tree->getOrder()+1,D))+7)/8);
    //NB: Gen nodes should take less space?
    this->sizeGenNodeCoeff = this->sizeNodeCoeff;
    println(10, "SizeNode Coeff (kB) " << this->sizeNodeCoeff*sizeof(double)/1024);
    println(10, "SizeGenNode Coeff (kB) " << this->sizeGenNodeCoeff*sizeof(double)/1024);

    //The first part of the Tree is filled with metadata; reserved size:
    this->sizeTreeMeta = 16*((sizeof(FunctionTree<D>)+15)/16);//we want only multiples of 16

    //FOR NOW WE DO NOT USE TREE METADATA
    this->sizeTreeMeta = 0;

    //The dynamical part of the tree is filled with nodes of size:
    this->sizeNode = this->sizeNodeMeta + this->sizeNodeCoeff*sizeof(double);//in units of Bytes
    this->sizeGenNode = this->sizeGenNodeMeta + this->sizeGenNodeCoeff*sizeof(double);//in units of Bytes

    //Tree is defined as array of doubles, because C++ does not like void malloc
    //NB: important to divide by sizeof(double) BEFORE multiplying. Otherwise we can get int overflow!
    int mysize = (this->sizeTreeMeta/sizeof(double) + this->maxNodes*(this->sizeNode/sizeof(double)));
    println(1, "Allocating array of size (MB)  " << mysize*sizeof(double)/1024/1024);
    this->SData = new double[mysize];
    //indicate occupation of nodes
    this->NodeStackStatus = new int[this->maxNodes+1];

    this->SGenData = new double[(this->maxGenNodes*sizeGenNode)/sizeof(double)];
    //indicate occupation of Gen nodes
    this->GenNodeStackStatus = new int[this->maxGenNodes+1];

    //NodeCoeff pointers
    this->CoeffStack = new double * [this->maxNodesCoeff];
    //indicate occupied NodeCoeff
    this->CoeffStackStatus = new int[this->maxNodesCoeff+1];
    //GenNodeCoeff pointers
    this->GenCoeffStack = new double * [this->maxGenNodesCoeff];
    //indicate occupied GenNodeCoeff
    this->GenCoeffStackStatus = new int[this->maxGenNodesCoeff+1];

    //some useful pointers to positions in SData
    this->lastNode = (ProjectedNode<D>*) (this->SData + this->sizeTreeMeta/sizeof(double));
    this->firstNode = (double*)this->lastNode;//constant start of node data
    this->lastNodeCoeff = (double*) (this->SData+this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNodeMeta/sizeof(double));//start after the metadata
    this->firstNodeCoeff = this->lastNodeCoeff;//constant start of coeff data

    this->lastGenNode = (GenNode<D>*) (this->SGenData);
    this->lastGenNodeCoeff = this->SGenData + (this->maxGenNodes*sizeGenNodeMeta)/sizeof(double);//start at after GenNode Meta part

    this->nNodesCoeff=-1;//add 1 before each allocation
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

    for (int i = 0; i <this->maxNodesCoeff;i++){
      this->CoeffStack[i] = this->lastNodeCoeff+i*this->sizeNodeCoeff;
      this->CoeffStackStatus[i] = 0;//0=unoccupied
    }
    this->CoeffStackStatus[this->maxNodesCoeff]=-1;//-1=unavailable

    for (int i = 0; i <this->maxGenNodesCoeff;i++){
      this->GenCoeffStack[i] = this->lastGenNodeCoeff+i*this->sizeGenNodeCoeff;
      this->GenCoeffStackStatus[i] = 0;//0=unoccupied
    }
    this->GenCoeffStackStatus[this->maxGenNodesCoeff] = -1;//-1=unavailable

    
    this->mwTree_p = Tree;//pointer to parent tree
    this->mwTree_p->serialTree_p = this;//must be defined before nodes are created

#ifdef HAVE_OPENMP
    omp_init_lock(&Stree_lock);
#endif

    //alloc root nodes
    MWNode<D> **roots = Tree->getRootBox().getNodes();
    int NodeIx;
    for (int rIdx = 0; rIdx < Tree->getRootBox().size(); rIdx++) {
        const NodeIndex<D> &nIdx = Tree->getRootBox().getNodeIndex(rIdx);
	roots[rIdx] = new (allocNodes(1, &NodeIx)) ProjectedNode<D>(*this->getTree(), nIdx);
	roots[rIdx]->NodeRank = NodeIx;
    }
    println(10, "tree at " << Tree<<" first root at "<<(Tree->getRootBox().getNodes())[0]<<" rootrank "<<(Tree->getRootBox().getNodes())[0]->getRank()<<" STree at "<<this);
    //cout<<MPI_rank<<" root 0 cvptr " << (double*)(*(char**)(roots[0]))<<endl;

    Tree->resetEndNodeTable();

}
/** Overwrite all pointers defined in the tree.
  * Necessary after sending the tree 
  * could be optimized. Should reset other counters? (GenNodes...) */
template<int D>
void SerialTree<D>::RewritePointers(){
  int DepthMax = 100;
  MWNode<D>* stack[DepthMax*8];
  int  slen = 0, counter = 0;

  //  ptrdiff_t d_p_shift = this->SData - *((double**)this->SData);//actual - written adress
  ptrdiff_t d_p_shift = 0;//actual - written adress of double
  ptrdiff_t n_p_shift = 0;//actual - written adress of node
  ptrdiff_t t_p_shift = 0;//actual - written adress of tree
  cout<<"start of data at "<<this->firstNode<<" pointer at "<<" diff "<<d_p_shift<<endl;
  cout<<" number of gen coeff on stack "<< this->nGenNodesCoeff <<" tag "<<this<<endl;
  cout<<" size tree "<< this->getTree()->nNodes <<" tag "<<this<<endl;

  this->lastNodeCoeff = (double*) (this->SData+this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNodeMeta/sizeof(double));//start after the metadata
  for (int i = 0; i <= this->nNodesCoeff; i++){
    CoeffStack[i] += d_p_shift;
  }

  FunctionTree<D>* Tree = dynamic_cast<FunctionTree<D>*> (this->mwTree_p);

  NodeBox<D> &rBox = Tree->getRootBox();
  MWNode<D> **roots = rBox.getNodes();
  //const NodeIndex<D> &nIdx = Tree->getRootBox().getNodeIndex(0);
  NodeIndex<D> nIdx;
  //ProjectedNode<D>* tempnode = new (allocNodes(1)) ProjectedNode<D>(*Tree, nIdx);//delete does not work for this
  ProjectedNode<D>* tempnode = new ProjectedNode<D>(*Tree, nIdx);
  //virtual table pointer. Must be made from a locally created node.
  char* cvptr =  *(char**)(tempnode);
  cout<<MPI_rank<<" Node virtual table pointer "<<(double*) cvptr <<endl;  
  delete tempnode; 
 
  d_p_shift = this->firstNodeCoeff-(*((double**) ((void*)&(roots[0]->coefvec))));
  n_p_shift =  ((MWNode<D> *) ( this->firstNodeCoeff)) - ((MWNode<D> *)(*((double**) ((void*)&(roots[0]->coefvec)))));
  //n_p_shift =  (d_p_shift*sizeof(double))/sizeof(MWNode<D>);//be careful with ptrdiff_t arithmetic!
  t_p_shift =  ((MWTree<D> *) ( this->firstNodeCoeff)) - ((MWTree<D> *)(*((double**) ((void*)&(roots[0]->coefvec)))));
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
      *(char**)(fpos) = cvptr;

      fpos->parent = (MWNode<D>*)(((double*)fpos->parent) + d_p_shift);// += n_p_shift not safe!
      //fpos->tree += t_p_shift;//NB: does not work!//(MWTree<D>*)(((double*)fpos->tree) + d_p_shift);
      fpos->tree = this->getTree();
      *((double**) ((void*)&(fpos->coefvec))) = this->firstNodeCoeff+(fpos->NodeCoeffIx)*this->sizeNodeCoeff;
      //*((int*) (((void*)&(fpos->coefvec))+8)) = this->tree->allocator->sizeNodeCoeff;

     fpos->coefs = &fpos->coefvec;
  }
  this->getTree()->resetEndNodeTable();

}
/** Adds two trees.
  * Generate missing nodes on the flight and adds all nodes */
template<int D>
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
    if(MPI_size == 2)SendRcv_SerialTree(TreeA, 0, 1, 44, MPI_COMM_WORLD);
#endif


}
/** Adds two trees.
  * Generate missing nodes on the flight and "compress" the ancestor  from summed nodes also on the flight */
template<int D>
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

}

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

  cA=&(Node->getCoefs()(0));
  cA=&(Node->children[0]->getCoefs()(0));

  double* coeffin  = &(Node->getCoefs()(0));
  double* coeffout = &(Node->children[0]->getCoefs()(0));

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
template<int D>
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

    
}

/** Make parent from children scaling coefficients
 * Other node info are not used/set
 * coeff_in is not modified.
 * The output is read directly from the 8 children scaling coefficients. 
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
template<int D>
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
      /* Operate in direction i only if the bits along other
       * directions are identical. The bit of the direction we
       * operate on determines the appropriate filter/operator */
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
  //  if(D%2)for(int i=0; i<kp1_d*tDim; i++) coeff_out[i] = tmpcoeff[i];

}

//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* SerialTree<D>::allocNodes(int nAlloc, int* NodeIx) {
  *NodeIx = this->nNodes;
  if (this->nNodes+nAlloc > this->maxNodes){
    println(0, "maxNodes exceeded " << this->maxNodes);
    MSG_FATAL("maxNodes exceeded ");
    this->nNodes -= nAlloc;
    return 0;
  } else {
      ProjectedNode<D>* newNode  = this->lastNode;
      ProjectedNode<D>* newNode_cp  = newNode;
      //We can have sizeNodeMeta with different size than ProjectedNode
      //println(0, "new size meta " << this->nNodes);
      for (int i = 0; i<nAlloc; i++){
	// DOES NOT WORK FOR LOOSE NODES: newNode_cp->NodeRank = *NodeIx+i;//Until overwritten!
	if (this->NodeStackStatus[*NodeIx+i]!=0)
	  println(0, *NodeIx+i<<" NodeStackStatus: not available " << this->NodeStackStatus[*NodeIx+i]);
	this->NodeStackStatus[*NodeIx+i]=1;
	newNode_cp++;
      }
      this->nNodes += nAlloc;
      this->lastNode = (ProjectedNode<D>*) ((char *) this->lastNode + nAlloc*this->sizeNodeMeta);
      return newNode;
  }
}
template<int D>
void SerialTree<D>::DeAllocNodes(int NodeRank) {
  if (this->nNodes <0){
    println(0, "minNodes exceeded " << this->nNodes);
    this->nNodes++;
  }
  this->NodeStackStatus[NodeRank]=0;//mark as available
  if(NodeRank==this->nNodes-1){//top of stack
    int TopStack=this->nNodes;
    while(this->NodeStackStatus[TopStack-1]==0){
      TopStack--;
      if(TopStack<1)break;
    }
    this->nNodes=TopStack;//move top of stack
  }
 }

//return pointer to the last active node or NULL if failed
template<int D>
GenNode<D>* SerialTree<D>::allocGenNodes(int nAlloc, int* NodeIx) {

  omp_set_lock(&Stree_lock);
  *NodeIx = this->nGenNodes;
  if (this->nGenNodes+nAlloc > this->maxGenNodes){
    println(0, "maxNodes exceeded " << this->maxGenNodes);
    MSG_FATAL("maxNodes exceeded ");
    omp_unset_lock(&Stree_lock);
    return 0;
  } else {
    GenNode<D>* newNode  = (GenNode<D>*)this->lastGenNode;
    GenNode<D>* newNode_cp  = newNode;
    for (int i = 0; i<nAlloc; i++){
      newNode_cp->NodeRank = *NodeIx+i;//Until overwritten!
      if (this->GenNodeStackStatus[*NodeIx+i]!=0)
	println(0, *NodeIx+i<<" NodeStackStatus: not available " << this->GenNodeStackStatus[*NodeIx+i]);
      this->GenNodeStackStatus[*NodeIx+i]=1;
      newNode_cp++;
   }
    this->nGenNodes += nAlloc;
    this->lastGenNode = (GenNode<D>*) ((char *) this->lastGenNode + nAlloc*this->sizeGenNodeMeta);
    //println(0, "new size meta " << this->nNodes);
    omp_unset_lock(&Stree_lock);
    return newNode;
  }
}

template<int D>
void SerialTree<D>::DeAllocGenNodes(int NodeRank) {
  omp_set_lock(&Stree_lock);
  if (this->nGenNodes <0){
    println(0, "minNodes exceeded " << this->nGenNodes);
    this->nGenNodes++;
  }
  this->GenNodeStackStatus[NodeRank]=0;//mark as available
  if(NodeRank==this->nGenNodes-1){//top of stack
    int TopStack=this->nGenNodes;
    while(this->GenNodeStackStatus[TopStack-1]==0){
      TopStack--;
      if(TopStack<1)break;
    }
    this->nGenNodes=TopStack;//move top of stack
  }
  omp_unset_lock(&Stree_lock);
 }

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocCoeff(int nAllocCoeff, MWNode<D>* node) {
//#pragma omp critical
  omp_set_lock(&Stree_lock);
  this->nNodesCoeff += 1;//nAllocCoeff;
  if (nAllocCoeff!=1<<D) MSG_FATAL("Only 2**D implemented now!");
  if (this->nNodesCoeff > this->maxNodesCoeff ){
    println(0, "maxNodesCoeff exceeded " << this->maxNodesCoeff);
    MSG_FATAL("maxNodesCoeff exceeded ");
    omp_unset_lock(&Stree_lock);
    return 0;
  } else if( CoeffStackStatus[this->nNodesCoeff]!=0){
    println(0, this->nNodesCoeff<<" CoeffStackStatus: not available " <<CoeffStackStatus[this->nNodesCoeff] );      
    omp_unset_lock(&Stree_lock);
    return 0;
  }else{
    double* NodeCoeff=this->CoeffStack[this->nNodesCoeff];
    this->CoeffStackStatus[this->nNodesCoeff]=1;
    node->NodeCoeffIx = this->nNodesCoeff;
    //node->NodeCoeffIx = node->NodeRank;
    omp_unset_lock(&Stree_lock);
    return NodeCoeff;
  }
}

//"Deallocate" the stack
template<int D>
void SerialTree<D>::DeAllocCoeff(int DeallocIx) {  
  omp_set_lock(&Stree_lock);
  if (this->CoeffStackStatus[DeallocIx]==0){
    println(0, "deleting already unallocated coeff " << DeallocIx);
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

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* SerialTree<D>::allocGenCoeff(int nAllocCoeff, MWNode<D>* node) {
  omp_set_lock(&Stree_lock);
  if (nAllocCoeff!=1<<D) MSG_FATAL("Only 2**D implemented now!");
  this->nGenNodesCoeff += 1;//nAllocCoeff;
  if (this->nGenNodesCoeff > this->maxGenNodesCoeff ){
    println(0, "maxGenNodesCoeff exceeded " << this->maxGenNodesCoeff);
    MSG_FATAL("maxGenNodesCoeff exceeded ");
    this->nGenNodesCoeff -= 1;//nAllocCoeff;
    omp_unset_lock(&Stree_lock);
    return 0;
    //  } else if( GenCoeffStackStatus[this->nGenNodesCoeff]!=0){
  } else if( GenCoeffStackStatus[node->NodeRank]!=0){
    println(0, node->NodeRank<<" "<<this->nGenNodesCoeff<<" GenCoeffStackStatus: not available " <<GenCoeffStackStatus[this->nGenNodesCoeff] );      
    omp_unset_lock(&Stree_lock);
    return 0;
  }else{
    //    double* GenNodeCoeff=this->GenCoeffStack[this->nGenNodesCoeff];
    //this->GenCoeffStackStatus[this->nGenNodesCoeff]=1;
    //node->NodeCoeffIx = this->nGenNodesCoeff;
    double* GenNodeCoeff=this->GenCoeffStack[node->NodeRank];
    this->GenCoeffStackStatus[node->NodeRank]=1;
    node->NodeCoeffIx = node->NodeRank;
    //for (int i = 0; i < this->sizeGenNodeCoeff; i++)coeffVectorXd[i]=0.0;
    omp_unset_lock(&Stree_lock);
    return GenNodeCoeff;
  }
}

//"Deallocate" the stack
template<int D>
void SerialTree<D>::DeAllocGenCoeff(int DeallocIx) {
  omp_set_lock(&Stree_lock);
  if (this->GenCoeffStackStatus[DeallocIx]==0){
    println(0, "deleting already unallocated Gencoeff " << DeallocIx);
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
    println(10, "~SerialTree");
    MWNode<D> **roots = this->getTree()->getRootBox().getNodes();
    for (int i = 0; i < this->getTree()->getRootBox().size(); i++) {
        ProjectedNode<D> *node = static_cast<ProjectedNode<D> *>(roots[i]);
        node->~ProjectedNode();
        roots[i] = 0;
    }
    //this->getTree()->MWTreeDestruct();
    delete[] this->SData;
    delete[] this->NodeStackStatus;

    delete[] this->SGenData;
    delete[] this->GenNodeStackStatus;

    delete[] this->CoeffStack;
    delete[] this->CoeffStackStatus;

    delete[] this->GenCoeffStack;
    delete[] this->GenCoeffStackStatus;
}

template class SerialTree<1>;
template class SerialTree<2>;
template class SerialTree<3>;
