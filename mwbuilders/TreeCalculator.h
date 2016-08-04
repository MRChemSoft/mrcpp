#ifndef TREECALCULATOR_H
#define TREECALCULATOR_H

#include "MWNode.h"
#include "mrcpp_declarations.h"

template<int D>
class TreeCalculator {
public:
    TreeCalculator() { }
    virtual ~TreeCalculator() { }

    virtual MWNodeVector* getInitialWorkVector(MWTree<D> &tree) const {
        return tree.copyEndNodeTable();
    }

    virtual void calcNodeVector(MWNodeVector &nodeVec) {
#pragma omp parallel shared(nodeVec)
{
        int nNodes = nodeVec.size();
	
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = *nodeVec[n];
            calcNode(node);
        }
}
    }
protected:
    virtual void calcNode(MWNode<D> &node) = 0;
};

#endif // TREECALCULATOR_H
