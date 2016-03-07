#ifndef MWADDER_H
#define MWADDER_H

#include <vector>

#include "TreeBuilder.h"

template<int D>
class MWAdder : public TreeBuilder<D> {
public:
    MWAdder(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        NOT_IMPLEMENTED_ABORT;
    }
    MWAdder(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        NOT_IMPLEMENTED_ABORT;
    }
    virtual ~MWAdder() {
        NOT_IMPLEMENTED_ABORT;
    }

    FunctionTree<D>* operator()(int a, FunctionTree<D> &a_tree,
                                int b, FunctionTree<D> &b_tree) {
        if (this->adaptor == 0) NOT_IMPLEMENTED_ABORT;
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, a, a_tree, b, b_tree);
        return out;
    }

    FunctionTree<D>* operator()(std::vector<int> coefs,
                                std::vector<FunctionTree<D> *> trees) {
        if (this->adaptor == 0) NOT_IMPLEMENTED_ABORT;
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, coefs, trees);
        return out;
    }

    void operator()(FunctionTree<D> &out,
                    int a, FunctionTree<D> &a_tree,
                    int b, FunctionTree<D> &b_tree) {
        NOT_IMPLEMENTED_ABORT;
    }

    void operator()(FunctionTree<D> &out,
                    std::vector<int> coefs,
                    std::vector<FunctionTree<D> *> trees) {
        NOT_IMPLEMENTED_ABORT;
    }
};

#endif // MWADDER_H
