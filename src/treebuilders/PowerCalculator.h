#pragma once

#include "TreeCalculator.h"

namespace mrcpp {

template<int D>
class PowerCalculator final : public TreeCalculator<D> {
public:
    PowerCalculator(FunctionTree<D> &inp, double pow) : power(pow), func(&inp) { }

private:
    double power;
    FunctionTree<D> *func;

    void calcNode(MWNode<D> &node_o) {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        int n_coefs = node_o.getNCoefs();
        double *coefs_o = node_o.getCoefs();
        // This generates missing nodes
        MWNode<D> node_i = func->getNode(idx); // Copy node
        node_i.mwTransform(Reconstruction);
        node_i.cvTransform(Forward);
        const double *coefs_i = node_i.getCoefs();
        for (int j = 0; j < n_coefs; j++) {
            coefs_o[j] = std::pow(coefs_i[j], this->power);
        }
        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

}
