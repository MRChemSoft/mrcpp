#ifndef MWMULTIPLIER_H
#define MWMULTIPLIER_H

#include <vector>

#include "TreeBuilder.h"

template<int D>
class MWMultiplier : public TreeBuilder<D> {
public:
    MWMultiplier(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        NOT_IMPLEMENTED_ABORT;
    }
    MWMultiplier(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        NOT_IMPLEMENTED_ABORT;
    }
    virtual ~MWMultiplier() {
        NOT_IMPLEMENTED_ABORT;
    }

    FunctionTree<D>* operator()(int coef,
                                FunctionTree<D> &a_tree,
                                FunctionTree<D> &b_tree) {
        if (this->adaptor == 0) NOT_IMPLEMENTED_ABORT;
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, coef, a_tree, b_tree);
        return out;
    }

    FunctionTree<D>* operator()(int coef, std::vector<FunctionTree<D> *> trees) {
        if (this->adaptor == 0) NOT_IMPLEMENTED_ABORT;
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, coef, trees);
        return out;
    }

    void operator()(FunctionTree<D> &out, int coef,
                    FunctionTree<D> &a_tree,
                    FunctionTree<D> &b_tree) {
        NOT_IMPLEMENTED_ABORT;
    }

    void operator()(FunctionTree<D> &out, int coef,
                    std::vector<FunctionTree<D> *> trees) {
        NOT_IMPLEMENTED_ABORT;
    }
};


#endif // MWMULTIPLIER_H
