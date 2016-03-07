#ifndef MRTREES_H_
#define MRTREES_H_

#include "MWTree.h"
//#include "FunctionTree.h"
#include "ProjectedNode.h"


/** Serialized MWTree.
 *  Contains all data to define tree and nodes, stored in a single contiguous piece of memory */

/** constructor.
 *  MaxSize is the max number of nodes the tree can have */
template<int D>
class MWTree_S{
 public:
  MWTree_S(const MultiResolutionAnalysis<D> &mra  ,int MaxNumberofNodes = 1 );
 
  ~MWTree_S();

  //The first part of the Tree is filled with metadata; reserved size:
  int SizeTreeMeta;
  //The dynamical part of the tree is filled with nodes of size:
  int SizeNode;

  //Tree is defined as array of doubles, because C++ does not like void malloc
  double* MWTree_S_array;

  MWTree<D>* MWTree_p;
  //  FunctionTree<D>* FunctionTree_p;

  int MaxNnodes;//max number of nodes that can be defined
  int Nnodes;//number of nodes already defined
  
  ProjectedNode<D>* LastNode;//pointer to the last active node

  ProjectedNode<D>* AllocNodes(int Nalloc);
  
  // Static default parameters
  const static int tDim = (1 << D);
};

#endif /* MRTREES_H_ */
