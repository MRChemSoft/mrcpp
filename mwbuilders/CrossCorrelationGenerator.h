#ifndef CROSSCORRELATIONGENERATOR_H
#define CROSSCORRELATIONGENERATOR_H

#include "TreeBuilder.h"
#include "OperatorAdaptor.h"
#include "CrossCorrelationCalculator.h"
#include "OperatorTree.h"
#include "Timer.h"

class CrossCorrelationGenerator : public TreeBuilder<2> {
public:
    CrossCorrelationGenerator(const MultiResolutionAnalysis<2> &mra, double pr)
            : TreeBuilder<2>(mra),
              prec(pr) {
    }
    virtual ~CrossCorrelationGenerator() {
    }

    void setPrecision(double pr) { this->prec = pr; }
    void multPrecision(double fac) { this->prec *= fac; }

    OperatorTree *operator()(FunctionTree<1> &inp) {
        OperatorTree *out = new OperatorTree(this->MRA, this->prec, MaxAllocNodes);
        (*this)(*out, inp, -1);
        return out;
    }

    void operator()(OperatorTree &out, FunctionTree<1> &inp, int maxIter = -1) {
        this->adaptor = new OperatorAdaptor(this->prec, this->MRA.getMaxScale());
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
