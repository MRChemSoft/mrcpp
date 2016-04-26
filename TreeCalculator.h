#ifndef TREECALCULATOR_H
#define TREECALCULATOR_H

#include "mwrepr_declarations.h"
#include "MWNode.h"

template<int D>
class TreeCalculator {
public:
    TreeCalculator() { }
    virtual ~TreeCalculator() { }

    virtual MWNodeVector* getInitialWorkVector(MWTree<D> &tree) const {
        return tree.copyEndNodeTable();
    }

    virtual double calcNodeVector(MWNodeVector &nodeVec) {
        double norm = 0.0;
        int nNodes = nodeVec.size();
#pragma omp parallel shared(nodeVec) firstprivate(nNodes) reduction(+:norm)
{
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = *nodeVec[n];
            calcNode(node);
            norm += node.getSquareNorm();
        }
}
        return std::max(norm, -1.0);
    }
protected:
    virtual void calcNode(MWNode<D> &node) = 0;
};

#endif // TREECALCULATOR_H
