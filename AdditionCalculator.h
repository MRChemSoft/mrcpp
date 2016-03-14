#ifndef ADDITIONCALCULATOR_H
#define ADDITIONCALCULATOR_H

#include "TreeCalculator.h"

template<int D>
class AdditionCalculator : public TreeCalculator<D> {
public:
    AdditionCalculator(std::vector<double> &c,
                       std::vector<FunctionTree<D> *> &f)
            : coefs(c),
              funcs(f) {
        if (c.size() != f.size()) MSG_ERROR("Invalid arguments");
    }
    virtual ~AdditionCalculator() { }

protected:
    std::vector<double> coefs;
    std::vector<FunctionTree<D> *> funcs;

    virtual void calcNode(MWNode<D> &outNode) const {
        const NodeIndex<D> &idx = outNode.getNodeIndex();
        Eigen::VectorXd &outVec = outNode.getCoefs();
        outNode.zeroCoefs();
        for (int n = 0; n < this->funcs.size(); n++) {
            const MWNode<D> &inpNode = this->funcs[n]->getNode(idx);
            const Eigen::VectorXd &inpVec = inpNode.getCoefs();
            outVec = outVec + this->coefs[n]*inpVec;
        }
        outNode.setHasCoefs();
        outNode.calcNorms();
    }
};

#endif // ADDITIONCALCULATOR_H
