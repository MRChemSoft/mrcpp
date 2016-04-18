#ifndef IDENTITYCALCULATOR_H
#define IDENTITYCALCULATOR_H

#include "TreeCalculator.h"
#include "MWNode.h"

template<int D>
class IdentityCalculator : public TreeCalculator<D> {
public:
    IdentityCalculator(FunctionTree<D> &inp_func) : func(&inp_func) { }
    virtual ~IdentityCalculator() { }

protected:
    FunctionTree<D> *func;
    virtual void calcNode(MWNode<D> &node) const {
        const NodeIndex<D> &idx = node.getNodeIndex();
        const MWNode<D> &inpNode = this->func->getNode(idx);
        node.setCoefs(inpNode.getCoefs());
        node.setHasCoefs();
        node.calcNorms();
    }
};

#endif // IDENTITYCALCULATOR_H
