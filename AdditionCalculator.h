#ifndef ADDITIONCALCULATOR_H
#define ADDITIONCALCULATOR_H

#include "TreeCalculator.h"

template<int D>
class AdditionCalculator : public TreeCalculator<D> {
public:
 AdditionCalculator(double aInp, FunctionTree<D> &a_funcInp,
                    double bInp, FunctionTree<D> &b_funcInp)
        : a(aInp),
          b(bInp),
          a_func(&a_funcInp),
          b_func(&b_funcInp) { }
    virtual ~AdditionCalculator() { }

protected:
    double a, b;
    FunctionTree<D> *a_func, *b_func;

    virtual void calcNode(MWNode<D> &node) const {
        const NodeIndex<D> &idx = node.getNodeIndex();
        const MWNode<D> &a_node = this->a_func->getNode(idx);
        const MWNode<D> &b_node = this->b_func->getNode(idx);
        const Eigen::VectorXd &aVec = a_node.getCoefs();
        const Eigen::VectorXd &bVec = b_node.getCoefs();
        Eigen::VectorXd &cVec = node.getCoefs();
        cVec = this->a * aVec + this->b * bVec;
        node.setHasCoefs();
        node.calcNorms();
     }
};



#endif // ADDITIONCALCULATOR_H
