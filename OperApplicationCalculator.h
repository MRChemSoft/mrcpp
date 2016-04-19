#ifndef OPERAPPLICATIONCALCULATOR_H
#define OPERAPPLICATIONCALCULATOR_H

#include "TreeCalculator.h"
#include "OperatorTreeVector.h"
#include "MWNode.h"

template<int D>
class OperApplicationCalculator : public TreeCalculator<D> {
public:
    OperApplicationCalculator(OperatorTreeVector &o, FunctionTree<D> &f)
        : oper(&o), func(&f) { }
    virtual ~OperApplicationCalculator() { }

protected:
    OperatorTreeVector *oper;
    FunctionTree<D> *func;

    virtual void calcNode(MWNode<D> &node) const {
        const NodeIndex<D> &idx = node.getNodeIndex();
        const MWNode<D> &inpNode = this->func->getNode(idx);
        const Eigen::VectorXd &inpVec = inpNode.getCoefs();
        Eigen::VectorXd &outVec = node.getCoefs();

        outVec = inpVec;
        node.setHasCoefs();
        node.calcNorms();
    }
};


#endif // OPERAPPLICATIONCALCULATOR_H
