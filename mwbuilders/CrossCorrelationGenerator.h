#ifndef CROSSCORRELATIONGENERATOR_H
#define CROSSCORRELATIONGENERATOR_H

#include "TreeBuilder.h"
#include "OperatorAdaptor.h"
#include "CrossCorrelationCalculator.h"
#include "OperatorTree.h"
#include "Timer.h"

class CrossCorrelationGenerator : public TreeBuilder<2> {
public:
    CrossCorrelationGenerator(double pr, int max_scale = MaxScale)
        : TreeBuilder<2>(pr, max_scale) { }
    virtual ~CrossCorrelationGenerator() { }

    void operator()(OperatorTree &out, FunctionTree<1> &inp, int maxIter = -1) {
        CrossCorrelationCalculator calculator(inp);
        OperatorAdaptor adaptor(this->prec, this->maxScale);
        this->build(out, calculator, adaptor, maxIter);

        Timer trans_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        out.setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
};

#endif // CROSSCORRELATIONGENERATOR_H
