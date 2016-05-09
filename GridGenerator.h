#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H

#include "TreeBuilder.h"
#include "DefaultCalculator.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"

template<int D>
class GridGenerator : public TreeBuilder<D> {
public:
    GridGenerator(const MultiResolutionAnalysis<D> &mra)
            : TreeBuilder<D>(mra) {
    }
    virtual ~GridGenerator() {
    }

    FunctionTree<D> *operator()() {
        return new FunctionTree<D>(this->MRA);
    }

    FunctionTree<D> *operator()(const RepresentableFunction<D> &inp,
                                int maxIter = -1) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp, maxIter);
        return out;
    }

    FunctionTree<D> *operator()(FunctionTree<D> &inp,
                                int maxIter = -1) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp, maxIter);
        return out;
    }

    FunctionTree<D> *operator()(FunctionTreeVector<D> &inp,
                                int maxIter = -1) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp, maxIter);
        return out;
    }

    void operator()(FunctionTree<D> &out,
                    const RepresentableFunction<D> &inp,
                    int maxIter = -1) {
        this->adaptor = new AnalyticAdaptor<D>(inp);
        this->calculator = new DefaultCalculator<D>();
        this->build(out, maxIter);
        this->clearCalculator();
        this->clearAdaptor();
        println(10, std::endl);
    }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int maxIter = -1) {
        this->adaptor = new CopyAdaptor<D>(inp);
        this->calculator = new DefaultCalculator<D>();
        this->build(out, maxIter);
        this->clearCalculator();
        this->clearAdaptor();
        println(10, std::endl);
    }

    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int maxIter = -1) {
        this->adaptor = new CopyAdaptor<D>(inp);
        this->calculator = new DefaultCalculator<D>();
        this->build(out, maxIter);
        this->clearCalculator();
        this->clearAdaptor();
        println(10, std::endl);
    }
};

#endif // GRIDGENERATOR_H
