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

#include <Eigen/Core>
#include "parallel.h"
#include "NodeIndex.h"
#include <vector>



template<int D> class MultiResolutionAnalysis;
template<int D> class ProjectedNode;
template<int D> class GenNode;
template<int D> class MWNode;
template<int D> class MWTree;
template<int D> class FunctionTree;
template<int D> class FunctionNode;

template<int D>
class SerialTree  {
public:
    SerialTree(MWTree<D>* Tree, int max_nodes);
    ~SerialTree();

    FunctionTree<D>* getTree() { return static_cast<FunctionTree<D> *>(this->mwTree_p); }

    ProjectedNode<D>* createSnode(const NodeIndex<D> &nIdx);
    ProjectedNode<D>* allocNodes(int Nalloc, int* NodeIx, double ** coefs_p);
    void DeAllocNodes(int NodeRank);
    GenNode<D>* allocGenNodes(int Nalloc, int* NodeIx, double ** coefs_p);
    void DeAllocGenNodes(int NodeRank);
    double* allocCoeff(int NallocCoeff, MWNode<D>* node);
    double* allocLooseCoeff(int NallocCoeff, MWNode<D>* node);
    void DeAllocLooseCoeff(int DeallocIx);
    double** CoeffStack;
    double** LooseCoeffStack;
    double** GenCoeffStack;
    void GenS_nodes(MWNode<D>* Node);
    void S_mwTransform(double* coeff_in, double* coeff_out, bool ReadOnlyScalingCoeff, int Children_Stride, bool overwrite=true);
    void S_mwTransformBack(double* coeff_in, double* coeff_out, int Children_Stride);

    void SerialTreeAdd(double c, FunctionTree<D>* &TreeB, FunctionTree<D>* &TreeC);
    void SerialTreeAdd_Up(double c, FunctionTree<D>* &TreeB, FunctionTree<D>* &TreeC);
    void RewritePointers(int Nchunks);
    int* NodeStackStatus;
    int* LooseNodeStackStatus;
    int* GenNodeStackStatus;
    int* LooseCoeffStackStatus;
    char* cvptr_ProjectedNode;//virtual table pointer for ProjectedNode
    char* cvptr_GenNode;// virtual table pointer for GenNode

    std::vector<ProjectedNode<D>*>  NodeChunks;
    std::vector<double*>  NodeCoeffChunks;
    std::vector<GenNode<D>*>  GenNodeChunks;
    std::vector<double*>  GenNodeCoeffChunks;
    int maxNodesPerChunk;

    friend class MWTree<D>;
    friend class ProjectedNode<D>;
    friend class MWNode<D>;
    friend class GenNode<D>;

    int nNodes;       //number of projected nodes already defined
    int nGenNodes;       //number of gen nodes already defined
    int nNodesCoeff;  //number of nodes Coeff already defined
    int nLooseNodesCoeff;  //number of loose nodes Coeff already defined
    int nGenNodesCoeff;  //number of Gen nodes Coeff already defined

    //    double* SData; //Nodes and coeff. Tree is defined as array of doubles, because C++ does not like void malloc
    ProjectedNode<D>* SNodes; //Serial Nodes 
    double* SNodesCoeff; //Serial Nodes coefficients
    double* LooseNodeCoeff; //To put coefficient of loose (temporary) nodes only
    //double* SGenData; //GenNodes and coeff
    GenNode<D>* SGenNodes;
    double* SGenNodesCoeff; //Serial Gen Nodes coefficients

    //    const MWTree<D>* mwTree_p;
    MWTree<D>* mwTree_p;
    ProjectedNode<D>* lastNode;//pointer to the last active node
    GenNode<D>* lastGenNode;//pointer to the last active Gen node
    int maxNodes;     //max number of nodes that can be defined
    int maxGenNodes;     //max number of Gen nodes that can be defined
    int maxLooseNodesCoeff;//max number of loose nodes Coeff that can be defined
    int sizeNodeCoeff;//size of coeff for one node
    int sizeGenNodeCoeff;//size of coeff for one Gen node
#ifdef HAVE_OPENMP
    omp_lock_t Stree_lock;
#endif
protected:
};

#endif /* TREEALLOCATOR_H_*/
