#ifndef OPERATORTREE_H
#define OPERATORTREE_H

#include "MWTree.h"

class OperatorTree: public MWTree<2> {
public:
    virtual ~OperatorTree();

protected:
    OperatorTree(const MultiResolutionAnalysis<2> &mra);
};

#endif // OPERATORTREE_H
