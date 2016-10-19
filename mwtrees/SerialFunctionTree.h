/**
*
*
*  \date Jul, 2016
*  \author Peter Wind <peter.wind@uit.no> \n
*  CTCC, University of Tromsø
*
*/

#ifndef SERIALFUNCTIONTREE_H_
#define SERIALFUNCTIONTREE_H_

#include <vector>

#include "SerialTree.h"
#include "parallel.h"

template<int D> class FunctionTree;
template<int D> class ProjectedNode;
template<int D> class GenNode;

template<int D>
class SerialFunctionTree : public SerialTree<D> {
public:
    SerialFunctionTree(FunctionTree<D> *tree, int max_nodes);
    virtual ~SerialFunctionTree();

    virtual void allocRoots(MWTree<D> &tree);
    virtual void allocChildren(MWNode<D> &parent);
    virtual void allocGenChildren(MWNode<D> &parent);

    virtual void deallocNodes(int serialIx);
    virtual void deallocGenNodes(int serialIx);

protected:
    int maxNodes;               //max number of nodes that can be defined
    int maxGenNodes;            //max number of Gen nodes that can be defined
    int sizeNodeCoeff;          //size of coeff for one node
    int sizeGenNodeCoeff;       //size of coeff for one Gen node

    int nNodes;                 //number of ProjectedNodes already defined
    int nGenNodes;              //number of GenNodes already defined
    int nNodesCoeff;            //number of nodes coeff already defined
    int nGenNodesCoeff;         //number of GenNodes coeff already defined

    double **coeffStack;
    double **genCoeffStack;

    int *nodeStackStatus;
    int *genNodeStackStatus;

    int maxNodesPerChunk;
    std::vector<ProjectedNode<D>*> nodeChunks;
    std::vector<GenNode<D>*> genNodeChunks;
    std::vector<double*> nodeCoeffChunks;
    std::vector<double*> genNodeCoeffChunks;

    ProjectedNode<D> *sNodes;   //serial ProjectedNodes
    GenNode<D> *sGenNodes;      //serial GenNodes
    double *sNodesCoeff;        //serial ProjectedNodes coefficients
    double *sGenNodesCoeff;     //serial GenNodes coefficients

    char *cvptr_ProjectedNode;  //virtual table pointer for ProjectedNode
    char *cvptr_GenNode;        //virtual table pointer for GenNode

    ProjectedNode<D>* lastNode; //pointer to the last active node
    GenNode<D>* lastGenNode;    //pointer to the last active Gen node

    ProjectedNode<D>* allocNodes(int nAlloc, int* serialIx, double **coefs_p);
    GenNode<D>* allocGenNodes(int nAlloc, int* serialIx, double **coefs_p);

private:
#ifdef HAVE_OPENMP
    omp_lock_t Sfunc_tree_lock;
#endif
};

#endif /* SERIALFUNCTIONTREE_H_*/
