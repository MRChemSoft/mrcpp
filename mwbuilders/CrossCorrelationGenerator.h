#ifndef CROSSCORRELATIONGENERATOR_H
#define CROSSCORRELATIONGENERATOR_H

#include "TreeBuilder.h"
#include "OperatorAdaptor.h"
#include "CrossCorrelationCalculator.h"
#include "OperatorTree.h"
#include "Timer.h"

class CrossCorrelationGenerator : public TreeBuilder<2> {
public:
    CrossCorrelationGenerator(double pr) : prec(pr) { }
    virtual ~CrossCorrelationGenerator() { }

    void setPrecision(double pr) { this->prec = pr; }
    void multPrecision(double fac) { this->prec *= fac; }

    void operator()(OperatorTree &out, FunctionTree<1> &inp, int maxIter = -1) {
        this->adaptor = new OperatorAdaptor(this->prec, MaxScale);
        this->calculator = new CrossCorrelationCalculator(inp);
        this->build(out, maxIter);
        this->clearCalculator();
        this->clearAdaptor();

        Timer trans_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        out.setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
protected:
    double prec;
};

#endif // CROSSCORRELATIONGENERATOR_H
