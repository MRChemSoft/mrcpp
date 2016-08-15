#include "TreeAllocator.h"
#include "MWNode.h"
#include "MWTree.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "ProjectedNode.h"
#include "GenNode.h"
#include "MathUtils.h"
#include "Timer.h"

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
    int SizeCoefOnly = (1<<D)*(MathUtils::ipow(mra.getOrder()+1,D))*sizeof(double);
    this->sizeNodeCoeff = SizeCoefOnly;
    //NB: Gen nodes should take less space?
    this->sizeGenNodeCoeff = SizeCoefOnly;
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
    println(0, "Allocating array of size (MB)  " << mysize*sizeof(double)/1024/1024);
    this->dataArray = new double[mysize];
    this->GenCoeffArray = new double[ this->maxGenNodesCoeff*this->sizeGenNodeCoeff/sizeof(double)];
    this->lastNode = (ProjectedNode<D>*) (this->dataArray + this->sizeTreeMeta/sizeof(double));
    //    this->lastNodeCoeff = (double*) (this->dataArray+this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNodeMeta/sizeof(double));
    this->lastNodeCoeff = (double*) (this->dataArray+this->sizeTreeMeta/sizeof(double) + this->maxNodes*this->sizeNodeMeta/sizeof(double));//start after the metadata
    this->firstNodeCoeff = this->lastNodeCoeff;//constant start of coeff data
    this->lastGenNodeCoeff = this->GenCoeffArray;//start at start of array
    this->nNodesCoeff=-1;//add 1 before each allocation
    this->nGenNodesCoeff=-1;//add 1 before each allocation

    NodeStackStatus = new int[maxNodes+1];
    for (int i = 0; i <maxNodes;i++){
      NodeStackStatus[i] = 0;//0=unoccupied
    }
    NodeStackStatus[maxNodes] = -1;//=unavailable

    CoeffStack = new double * [maxNodesCoeff];
    CoeffStackStatus = new int[maxNodesCoeff+1];
    for (int i = 0; i <maxNodesCoeff;i++){
      CoeffStack[i] = this->lastNodeCoeff+i*this->sizeNodeCoeff/sizeof(double);
      CoeffStackStatus[i] = 0;//0=unoccupied
    }
    println(0, "start of coeff stack at "<<CoeffStack[0]);
    CoeffStackStatus[maxNodesCoeff]=-1;//-1=unavailable
    //CoeffStackStatus[-1]=-1;//TO FIX!

    GenCoeffStack = new double * [maxGenNodesCoeff];
    GenCoeffStackStatus = new int[maxGenNodesCoeff+1];
    for (int i = 0; i <maxGenNodesCoeff;i++){
      GenCoeffStack[i] = this->lastGenNodeCoeff+i*this->sizeGenNodeCoeff/sizeof(double);
      GenCoeffStackStatus[i] = 0;//0=unoccupied
    }
    GenCoeffStackStatus[maxGenNodesCoeff] = -1;//-1=unavailable

     //put a MWTree at start of array
    this->mwTree_p = new (this->dataArray) MWTree<D>(mra);

    this->mwTree_p->allocator = this;

    //alloc root nodes
    NodeBox<D> &rBox = this->mwTree_p->getRootBox();
    MWNode<D> **roots = rBox.getNodes();
    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
        const NodeIndex<D> &nIdx = rBox.getNodeIndex(rIdx);
         roots[rIdx] = new (allocNodes(1)) ProjectedNode<D>(*this->getTree(), nIdx);
    }

    this->mwTree_p->resetEndNodeTable();
}

template<int D>
void TreeAllocator<D>::SerialTreeAdd(double c, FunctionTree<D>* &TreeB, FunctionTree<D>* &TreeC){
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
    int Children_Stride = this->sizeNodeCoeff/sizeof(double);
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
    timer.restart();
    t2.restart();
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
	  t1.restart();
	  GenS_nodes(fposA);
	  t1.stop();
	}
	if(fposB->getNChildren()==0){
	  //println(0, slenB<<" BgenChildren()  ");
	  t1.restart();
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
	      	      
	      t2.restart();
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
	      counterA++;
	    }
	  }
	  if(slenA>0){
	    //make parent
	    t3.restart();
	    S_mwTransformBack(&(fposA->getCoefs()(0)), &(fposA->parent->getCoefs()(0)), Children_Stride); 
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

    /*
   slenA=0;
   counterA=0;
   slenB=0;
   counterB=0;
    for (int rIdxA = 0; rIdxA < rBoxA.size(); rIdxA++) {
       const NodeIndex<D> &nIdx = rBoxA.getNodeIndex(rIdxA);
	stackA[slenA++] = TreeA->findNode(nIdx);
    }

    MWNode<D>* stackC[DepthMax*8];
    NodeBox<D> &rBoxC = TreeC->getRootBox();
    MWNode<D> **rootsC = rBoxC.getNodes();
    int slenC=0;
    int counterC=0;
    for (int rIdxC = 0; rIdxC < rBoxC.size(); rIdxC++) {
      const NodeIndex<D> &nIdx = rBoxC.getNodeIndex(rIdxC);
      stackC[slenC++] = TreeC->findNode(nIdx);
    }
    while (slenC) {
      counterC++;
      MWNode<D>* fposC = stackC[--slenC];
      MWNode<D>* fposA = stackA[--slenA];
     if(fposC->getNChildren()){
	for (int i = 0; i < fposC->getNChildren(); i++) {
	  stackC[slenC++] = fposC->children[i];
	  stackA[slenA++] = fposA->children[i];
	}
      }else{
	  cA=&(fposA->getCoefs()(0));
	  cC=&(fposC->getCoefs()(0));
	  if(fabs(cA[0]-cC[0])>4.E-3){
	    double* val=fposC->tree->allocator->firstNodeCoeff+1598*4096;
	    println(0, " indexA "<<fposA->getNodeIndex()<<" rankA "<<fposA->getRank()<<" indexC "<<fposC->getNodeIndex()<<" rankC "<<fposC->getRank()<<" coeffA "<<cA[0]<<" coeffC "<<cC[0]<<" at "<<&(cC[0])<<" coefix "<<fposC->NodeCoeffIx<<" val "<<val[0]);
	  }
      }
    }
    println(0, " Nodes counterC "<<counterC);

    */

    this->getTree()->squareNorm=Tsum;
    println(0, " time generate     " << t1);
    println(0, " time add coef     " << t2);
    t3.restart();
    //this->S_mwTreeTransformUp();
     //this->mwTree_p->mwTransform(BottomUp);
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
void TreeAllocator<D>::GenS_nodes(MWNode<D>* Node){

  bool ReadOnlyScalingCoeff=false;

  if(Node->isGenNode())ReadOnlyScalingCoeff=true;
  double* cA;

  Node->genChildren();//will make children and allocate coeffs, but without setting values for coeffs.

  cA=&(Node->getCoefs()(0));
  cA=&(Node->children[0]->getCoefs()(0));

  double* coeffin  = &(Node->getCoefs()(0));
  double* coeffout = &(Node->children[0]->getCoefs()(0));

  int Children_Stride = this->sizeGenNodeCoeff/sizeof(double);
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
void TreeAllocator<D>::S_mwTransform(double* coeff_in, double* coeff_out, bool ReadOnlyScalingCoeff, int Children_Stride) {
  if(D!=3)MSG_FATAL("S_mwtransform Only D=3 implemented now!");
  int operation = Reconstruction;
  int kp1 = this->mwTree_p->getKp1();
  int tDim = (1<<D);
  int kp1_dm1 = kp1*kp1;//MathUtils::ipow(Parent->getKp1(), D - 1)
  int kp1_d = kp1_dm1*kp1;//
  const MWFilter &filter = this->mwTree_p->getMRA().getFilter();
  //VectorXd &result = Parent->tree->getTmpMWCoefs();
  double overwrite = 0.0;
  double *tmp;
  double tmpcoeff[kp1_d*tDim];
  //double tmp2coeff[kp1_d*tDim];
  int ftlim=tDim;
  int ftlim2=tDim;
  int ftlim3=tDim;
  if(ReadOnlyScalingCoeff){
    ftlim=1;
    ftlim2=2;
    ftlim3=4;
    //NB: Careful: tmpcoeff and tmp2coeff are not initialized to zero
    //must not read these unitialized values either!
  }

  int i = 0;
  int mask = 1;
  for (int gt = 0; gt < tDim; gt++) {
    double *out = coeff_out + gt * kp1_d;
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
    //double *out = coeff_out + gt * kp1_d;
    double *out = coeff_out + gt * Children_Stride;
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

}


/** Regenerate all s/d-coeffs by backtransformation, starting at the bottom and
  * thus purifying all coefficients.
  * not yet fully optimized.
 */
template<int D>
void TreeAllocator<D>::S_mwTreeTransformUp() {
    Timer t0,t1,t2,t3;
    int DepthMax = 100;
    MWTree<D>* Tree=this->mwTree_p;
    int status[this->nNodes];
    MWNode<D>* stack[DepthMax*8];
    int  slen = 0, counter = 0, counter2 = 0;
    NodeBox<D> &rBox = Tree->getRootBox();
    MWNode<D> **roots = rBox.getNodes();
    int Children_Stride = this->sizeGenNodeCoeff/sizeof(double);
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
	    t0.restart();
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
void TreeAllocator<D>::S_mwTransformBack(double* coeff_in, double* coeff_out, int Children_Stride) {
  if(D!=3)MSG_FATAL("S_mwtransform Only D=3 implemented now!");
  int operation = Compression;
  int kp1 = this->mwTree_p->getKp1();
  int tDim = (1<<D);
  int kp1_dm1 = kp1*kp1;//MathUtils::ipow(Parent->getKp1(), D - 1)
  int kp1_d = kp1_dm1*kp1;//
  const MWFilter &filter = this->mwTree_p->getMRA().getFilter();
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
ProjectedNode<D>* TreeAllocator<D>::allocNodes(int nAlloc) {
  this->nNodes += nAlloc;
  if (this->nNodes > this->maxNodes){
    println(0, "maxNodes exceeded " << this->maxNodes);
    MSG_FATAL("maxNodes exceeded ");
    this->nNodes -= nAlloc;
    return 0;
  } else {
      ProjectedNode<D>* oldlastNode  = this->lastNode;
      //We can have sizeNodeMeta with different size than ProjectedNode
      this->lastNode = (ProjectedNode<D>*) ((char *) this->lastNode + nAlloc*this->sizeNodeMeta);
      //println(0, "new size meta " << this->nNodes);
      for (int i = 0; i<nAlloc; i++){
	oldlastNode->NodeRank=this->nNodes+i-1;
	if (this->NodeStackStatus[this->nNodes+i-1]!=0)
	  println(0, this->nNodes+i-1<<" NodeStackStatus: not available " << this->NodeStackStatus[this->nNodes+i-1]);
	this->NodeStackStatus[this->nNodes+i-1]=1;
      }
      return oldlastNode;
  }
}
template<int D>
void TreeAllocator<D>::DeAllocNodes(int NodeRank) {
  if (this->nNodes <0){
    println(0, "minNodes exceeded " << this->nNodes);
    this->nNodes++;
  }
  this->NodeStackStatus[NodeRank]=0;//mark as available
  if(NodeRank==this->nNodes-1){//top of stack
    int TopStack=this->nNodes;
    while(this->NodeStackStatus[TopStack-1]==0 and TopStack>0){
      TopStack--;
    }
    this->nNodes=TopStack;//move top of stack
  }
 }
//return pointer to the last active node or NULL if failed
template<int D>
GenNode<D>* TreeAllocator<D>::allocGenNodes(int nAlloc) {
    this->nNodes += nAlloc;
    if (this->nNodes > this->maxNodes){
      println(0, "maxNodes exceeded " << this->maxNodes);
      MSG_FATAL("maxNodes exceeded ");
       this->nNodes -= nAlloc;
        return 0;
    } else {
      GenNode<D>* oldlastNode  = (GenNode<D>*)this->lastNode;
	this->lastNode = (ProjectedNode<D>*) ((char *) this->lastNode + nAlloc*this->sizeNodeMeta);
	//We can have sizeNodeMeta with different size than ProjectedNode
	for (int i = 0; i<nAlloc; i++){
	  oldlastNode->NodeRank=this->nNodes+i-1;
	  if (this->NodeStackStatus[this->nNodes+i-1]!=0)
	    println(0, this->nNodes+i-1<<" NodeStackStatus: not available " << this->NodeStackStatus[this->nNodes+i-1]);
	  this->NodeStackStatus[this->nNodes+i-1]=1;
	}
        //println(0, "new size meta " << this->nNodes);
        return oldlastNode;
    }
}
//return pointer to the Coefficients of the node or NULL if failed
template<int D>
double* TreeAllocator<D>::allocCoeff(int nAllocCoeff) {
//VectorXd* TreeAllocator<D>::allocCoeff(int nAllocCoeff) {
  //for now only scaling and wavelets, because we need to take into account the size of vector class
  this->nNodesCoeff += 1;//nAllocCoeff;
    if (nAllocCoeff!=8) MSG_FATAL("Only 2**D implemented now!");
    if (this->nNodesCoeff > this->maxNodesCoeff ){
      println(0, "maxNodesCoeff exceeded " << this->maxNodesCoeff);
      MSG_FATAL("maxNodesCoeff exceeded ");
        return 0;
    } else if( CoeffStackStatus[this->nNodesCoeff]!=0){
      println(0, this->nNodesCoeff<<" CoeffStackStatus: not available " <<CoeffStackStatus[this->nNodesCoeff] );      
        return 0;
    }else{
      double* coeffVectorXd=this->CoeffStack[this->nNodesCoeff];
      this->CoeffStackStatus[this->nNodesCoeff]=1;
      for (int i = 0; i < this->sizeNodeCoeff/8; i++)coeffVectorXd[i]=0.0;
      //println(0,this->nNodesCoeff<< " Coeff returned: " << (double*)coeffVectorXd);      
      return coeffVectorXd;
    }
}

//"Deallocate" the stack
template<int D>
void TreeAllocator<D>::DeAllocCoeff(int DeallocIx) {
    if (this->CoeffStackStatus[DeallocIx]==0)println(0, "deleting already unallocated coeff " << DeallocIx);
    this->CoeffStackStatus[DeallocIx]=0;//mark as available
    if(DeallocIx==this->nNodesCoeff){//top of stack
      int TopStack=this->nNodesCoeff;
      while(this->CoeffStackStatus[TopStack]==0 and TopStack>0){
	TopStack--;
      }
      this->nNodesCoeff=TopStack;//move top of stack
    }
}

//return pointer to the Coefficients of the node or NULL if failed
template<int D>
//VectorXd* TreeAllocator<D>::allocGenCoeff(int nAllocCoeff) {
double* TreeAllocator<D>::allocGenCoeff(int nAllocCoeff) {
  //for now only scaling and wavelets, because we need to take into account the size of vector class
    if (nAllocCoeff!=8) MSG_FATAL("Only 2**D implemented now!");
    this->nGenNodesCoeff += 1;//nAllocCoeff;
    if (this->nGenNodesCoeff > this->maxGenNodesCoeff ){
      println(0, "maxGenNodesCoeff exceeded " << this->maxGenNodesCoeff);
       MSG_FATAL("maxGenNodesCoeff exceeded ");
      this->nGenNodesCoeff -= 1;//nAllocCoeff;
        return 0;
    } else if( GenCoeffStackStatus[this->nGenNodesCoeff]!=0){
      println(0, this->nGenNodesCoeff<<" GenCoeffStackStatus: not available " <<GenCoeffStackStatus[this->nGenNodesCoeff] );      
        return 0;
    }else{
      double* coeffVectorXd=this->GenCoeffStack[this->nGenNodesCoeff];
      this->GenCoeffStackStatus[this->nGenNodesCoeff]=1;
      for (int i = 0; i < this->sizeGenNodeCoeff/8; i++)coeffVectorXd[i]=0.0;
      return coeffVectorXd;
    }
}

//"Deallocate" the stack
template<int D>
void TreeAllocator<D>::DeAllocGenCoeff(int DeallocIx) {
    if (this->GenCoeffStackStatus[DeallocIx]==0){
      println(0, "deleting already unallocated Gencoeff " << DeallocIx);
    }else{
      this->GenCoeffStackStatus[DeallocIx]=0;//mark as available
      if(DeallocIx==this->nGenNodesCoeff){//top of stack
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
    //    for (int i = 0; i <maxNodesCoeff;i++)this->CoeffStack[i]->~VectorXd();
    //    for (int i = 0; i <maxGenNodesCoeff;i++)this->GenCoeffStack[i]->~VectorXd();
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
