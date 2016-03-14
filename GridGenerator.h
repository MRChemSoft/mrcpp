#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H

#include "TreeBuilder.h"
#include "DefaultCalculator.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"

template<int D>
class GridGenerator : public TreeBuilder<D> {
public:
    GridGenerator(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->calculator = new DefaultCalculator<D>();
    }
    virtual ~GridGenerator() {
        this->clearCalculator();
    }

    FunctionTree<D> *operator()() {
        return new FunctionTree<D>(this->MRA);
    }

    FunctionTree<D> *operator()(const RepresentableFunction<D> &inp) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp);
        return out;
    }

    FunctionTree<D> *operator()(const FunctionTree<D> &inp) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp);
        return out;
    }

    void operator()(FunctionTree<D> &out, const RepresentableFunction<D> &inp) {
        this->adaptor = new AnalyticAdaptor<D>(inp);
        this->build(out);
        this->clearAdaptor();
        println(10, std::endl);
    }

    void operator()(FunctionTree<D> &out, FunctionTree<D> &inp) {
        this->adaptor = new CopyAdaptor<D>(inp);
        this->build(out);
        this->clearAdaptor();
        println(10, std::endl);
    }
};

#endif // GRIDGENERATOR_H
