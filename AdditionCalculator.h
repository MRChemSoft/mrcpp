#ifndef ADDITIONCALCULATOR_H
#define ADDITIONCALCULATOR_H

#include "TreeCalculator.h"
#include "FunctionTreeVector.h"

template<int D> class MWAdder;

template<int D>
class AdditionCalculator : public TreeCalculator<D> {
public:
    friend class MWAdder<D>;
protected:
    AdditionCalculator(FunctionTreeVector<D> &inp) : sum_vec(&inp) { }
    virtual ~AdditionCalculator() { }

    virtual void calcNode(MWNode<D> &node_o) const {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        Eigen::VectorXd &vec_o = node_o.getCoefs();
        vec_o.setZero();
        for (int i = 0; i < this->sum_vec->size(); i++) {
            double coef_i = this->sum_vec->getCoef(i);
            FunctionTree<D> &func_i = this->sum_vec->getFunc(i);
            const MWNode<D> &node_i = func_i.getNode(idx);
            const Eigen::VectorXd &vec_i = node_i.getCoefs();
            vec_o = vec_o + coef_i*vec_i;
        }
        node_o.setHasCoefs();
        node_o.calcNorms();
    }

private:
    FunctionTreeVector<D> *sum_vec;
};

#endif // ADDITIONCALCULATOR_H
