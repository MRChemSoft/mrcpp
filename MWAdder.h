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
        this->adaptor = a.copy();
    }
    virtual ~MWAdder() { this->clearAdaptor(); }

    FunctionTree<D>* operator()(int a, FunctionTree<D> &a_tree,
                                int b, FunctionTree<D> &b_tree) {
        std::vector<double> coefs;
        std::vector<FunctionTree<D> *> trees;
        coefs.push_back(a);
        coefs.push_back(b);
        trees.push_back(&a_tree);
        trees.push_back(&b_tree);
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, coefs, trees);
        return out;
    }

    FunctionTree<D>* operator()(std::vector<double> coefs,
                                std::vector<FunctionTree<D> *> trees) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, coefs, trees);
        return out;
    }

    void operator()(FunctionTree<D> &out,
                    int a, FunctionTree<D> &a_tree,
                    int b, FunctionTree<D> &b_tree) {
        std::vector<double> coefs;
        std::vector<FunctionTree<D> *> trees;
        coefs.push_back(a);
        coefs.push_back(b);
        trees.push_back(&a_tree);
        trees.push_back(&b_tree);
        (*this)(out, coefs, trees);
    }

    void operator()(FunctionTree<D> &out,
                    std::vector<double> coefs,
                    std::vector<FunctionTree<D> *> trees) {
        bool defaultAdaptor = false;
        if (this->adaptor == 0) {
            defaultAdaptor = true;
            this->adaptor = new CopyAdaptor<D>(trees);
        }
        this->calculator = new AdditionCalculator<D>(coefs, trees);
        this->build(out);
        out.mwTransform(BottomUp);
        for (int n = 0; n < trees.size(); n++) {
            trees[n]->deleteGenerated();
        }
        this->clearCalculator();
        if (defaultAdaptor) this->clearAdaptor();
    }
};

#endif // MWADDER_H
