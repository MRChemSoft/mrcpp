#pragma once

#include "TreeCalculator.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {

template<int D>
class AdditionCalculator final : public TreeCalculator<D> {
public:
    AdditionCalculator(const FunctionTreeVector<D> &inp) : sum_vec(inp) { }

protected:
    FunctionTreeVector<D> sum_vec;

    void calcNode(MWNode<D> &node_o) {
        node_o.zeroCoefs();
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        double *coefs_o = node_o.getCoefs();
        for (int i = 0; i < this->sum_vec.size(); i++) {
            double c_i = get_coef(this->sum_vec, i);
            FunctionTree<D> &func_i = get_func(this->sum_vec, i);
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
