#pragma once

#include "mwbuilders/TreeCalculator.h"
#include "mwtrees/FunctionTreeVector.h"

namespace mrcpp {

template<int D>
class AdditionCalculator : public TreeCalculator<D> {
public:
    AdditionCalculator(FunctionTreeVector<D> &inp) : sum_vec(&inp) { }
    virtual ~AdditionCalculator() { }

protected:
    FunctionTreeVector<D> *sum_vec;

    virtual void calcNode(MWNode<D> &node_o) {
        node_o.zeroCoefs();
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        double *coefs_o = node_o.getCoefs();
        for (int i = 0; i < this->sum_vec->size(); i++) {
            double c_i = this->sum_vec->getCoef(i);
            FunctionTree<D> &func_i = this->sum_vec->getFunc(i);
            // This generates missing nodes
            const MWNode<D> &node_i = func_i.getNode(idx);
            const double *coefs_i = node_i.getCoefs();
            int n_coefs = node_i.getNCoefs();
            for (int j = 0; j < n_coefs; j++) {
                coefs_o[j] += c_i * coefs_i[j];
            }
        }
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

}
