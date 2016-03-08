#ifndef MWADDER_H
#define MWADDER_H

#include <vector>

#include "TreeBuilder.h"
#include "AdditionCalculator.h"

template<int D>
class MWAdder : public TreeBuilder<D> {
public:
    MWAdder(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
    }
    MWAdder(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        NOT_IMPLEMENTED_ABORT;
    }
    virtual ~MWAdder() {
    }

    FunctionTree<D>* operator()(int a, FunctionTree<D> &a_tree,
                                int b, FunctionTree<D> &b_tree) {
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
        this->calculator = new AdditionCalculator<D>(a, a_tree, b, b_tree);
        this->build(out);
        out.mwTransform(BottomUp);
        this->clearCalculator();
    }

    void operator()(FunctionTree<D> &out,
                    std::vector<int> coefs,
                    std::vector<FunctionTree<D> *> trees) {
        NOT_IMPLEMENTED_ABORT;
    }
};

#endif // MWADDER_H
