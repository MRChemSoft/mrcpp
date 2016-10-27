#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H

#include "TreeBuilder.h"
#include "DefaultCalculator.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"
#include "SerialTree.h"

template<int D>
class GridGenerator : public TreeBuilder<D> {
public:
    GridGenerator() { }
    virtual ~GridGenerator() { }

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
