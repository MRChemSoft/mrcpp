#ifndef DEFAULTCALCULATOR_H
#define DEFAULTCALCULATOR_H

#include "TreeCalculator.h"

template<int D>
class DefaultCalculator :public TreeCalculator<D> {
public:
    DefaultCalculator() { }
    virtual ~DefaultCalculator() { }

    // Reimplementation without OpenMP, the default is faster this way
    virtual double calcNodeVector(MWNodeVector &nodeVec) const {
        int nNodes = nodeVec.size();
        for (int n = 0; n < nNodes; n++) {
            calcNode(*nodeVec[n]);
        }
        return -1.0;
    }
protected:
    virtual void calcNode(MWNode<D> &node) const {
        node.clearHasCoefs();
        node.clearNorms();
    }
};

#endif // DEFAULTCALCULATOR_H
