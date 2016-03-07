#include "MWTree.h"
#include "MWTree_S.h"

using namespace std;

/** Serialized MWTree.
 *  Contains all data to define tree and nodes, stored in a single contiguous piece of memory
 *  Maxnodes is the max number of nodes the tree can have */
template<int D>
MWTree_S<D>::MWTree_S(const MultiResolutionAnalysis<D> &mra , int MaxNumberofNodes) {
  
  //The first part of the Tree is filled with metadata; reserved size:
  SizeTreeMeta = (sizeof(MWTree<D>)+7)/sizeof(double);
  //The dynamical part of the tree is filled with nodes of size:
  SizeNode = (sizeof(ProjectedNode<D>)+7)/sizeof(double);
  cout<<"SizeNode "<<sizeof(ProjectedNode<D>)<<endl;

  //Tree is defined as array of doubles, because C++ does not like void malloc
  MWTree_S_array = new double[SizeTreeMeta + MaxNumberofNodes*SizeNode];

  MWTree_p = (MWTree<D>*) MWTree_S_array;
  //  FunctionTree_p = ( FunctionTree<D>*) MWTree_S_array;

  MaxNnodes = MaxNumberofNodes;//max number of nodes that can be defined
  Nnodes = 0;//number of nodes already defined
  
  LastNode = (ProjectedNode<D>*) (MWTree_S_array + SizeTreeMeta + Nnodes*SizeNode);

  MWTree<D> *  f_tree_S_start = new (MWTree_S_array) FunctionTree<D>(mra);//put a MWTree at start of array
  
  // Static default parameters
  const static int tDim = (1 << D);
}

//return pointer to the last active node or NULL if failed
template<int D>
ProjectedNode<D>* MWTree_S<D>::AllocNodes(int Nalloc) {    
  Nnodes += Nalloc;
  if(Nnodes > MaxNnodes){
    Nnodes -= Nalloc;
    return NULL;
  }else{    
    LastNode += Nalloc*SizeNode;
    cout<<"new size "<<Nnodes<<endl;
    return LastNode-Nalloc*SizeNode;
  }
}
/** MWTree_S destructor. */
template<int D>
MWTree_S<D>::~MWTree_S() {
  //  delete this->FunctionTree_p;
  delete this->MWTree_S_array;

}

template class MWTree_S<1>;
template class MWTree_S<2>;
template class MWTree_S<3>;
