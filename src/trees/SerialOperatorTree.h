/**
*
*
*  \date Jul, 2016
*  \author Peter Wind <peter.wind@uit.no> \n
*  CTCC, University of Tromsø
*
*/

#pragma once

#include <vector>

#include "SerialTree.h"

namespace mrcpp {

class SerialOperatorTree final : public SerialTree<2> {
public:
    SerialOperatorTree(OperatorTree *tree);
    SerialOperatorTree(const SerialOperatorTree &tree) = delete;
    SerialOperatorTree &operator=(const SerialOperatorTree &tree) = delete;
    ~SerialOperatorTree();

    void allocRoots(MWTree<2> &tree);
    void allocChildren(MWNode<2> &parent);
    void allocGenChildren(MWNode<2> &parent);

    void deallocNodes(int serialIx);
    void deallocGenNodes(int serialIx);
    void deallocGenNodeChunks();

protected:
    OperatorNode *sNodes;       //serial OperatorNodes

    std::vector<OperatorNode*> nodeChunks;
    std::vector<double*> nodeCoeffChunks;

    char *cvptr_OperatorNode;   //virtual table pointer for OperatorNode
    OperatorNode* lastNode;     //pointer to the last active node

    OperatorNode* allocNodes(int nAlloc, int* serialIx, double **coefs_p);

private:
#ifdef HAVE_OPENMP
    omp_lock_t Soper_tree_lock;
#endif
};

}
