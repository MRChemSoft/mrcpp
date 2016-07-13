/**
*
*
*  \date Jul, 2016
*  \author Peter Wind <peter.wind@uit.no> \n
*  CTCC, University of Tromsø
*
*/

#ifndef TREEALLOCATOR_H_
#define TREEALLOCATOR_H_

template<int D> class MultiResolutionAnalysis;
template<int D> class ProjectedNode;
template<int D> class MWTree;
template<int D> class FunctionTree;
template<int D> class FunctionNode;

template<int D>
class TreeAllocator  {
public:
    TreeAllocator(const MultiResolutionAnalysis<D> &mra, int max_nodes);
    virtual ~TreeAllocator();

    FunctionTree<D>* getTree() { return static_cast<FunctionTree<D> *>(this->mwTree_p); }

    ProjectedNode<D>* allocNodes(int Nalloc);

protected:
    int sizeTreeMeta; //The first part of the Tree is filled with metadata; reserved size:
    int sizeNode;     //The dynamical part of the tree is filled with nodes of size:
    int maxNodes;     //max number of nodes that can be defined
    int nNodes;       //number of nodes already defined

    double* dataArray; //Tree is defined as array of doubles, because C++ does not like void malloc

    MWTree<D>* mwTree_p;
    ProjectedNode<D>* lastNode;//pointer to the last active node
};

#endif /* TREEALLOCATOR_H_*/
