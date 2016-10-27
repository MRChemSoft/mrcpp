#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H

#include "TreeBuilder.h"
#include "DefaultCalculator.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"

template<int D>
class GridGenerator : public TreeBuilder<D> {
public:
    GridGenerator(int max_scale = MaxScale) : TreeBuilder<D>(-1.0, max_scale) { }
    virtual ~GridGenerator() { }

    void operator()(FunctionTree<D> &out,
                    const RepresentableFunction<D> &inp,
                    int maxIter = -1) {
        AnalyticAdaptor<D> adaptor(inp);
        DefaultCalculator<D> calculator;
        this->build(out, calculator, adaptor, maxIter);
        println(10, std::endl);
    }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int maxIter = -1) {
        CopyAdaptor<D> adaptor(inp);
        DefaultCalculator<D> calculator;
        this->build(out, calculator, adaptor, maxIter);
        println(10, std::endl);
    }

    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int maxIter = -1) {
        CopyAdaptor<D> adaptor(inp);
        DefaultCalculator<D> calculator;
        this->build(out, calculator, adaptor, maxIter);
        println(10, std::endl);
    }
};

#endif // GRIDGENERATOR_H
