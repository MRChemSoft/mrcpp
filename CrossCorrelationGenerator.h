#ifndef CROSSCORRELATIONGENERATOR_H
#define CROSSCORRELATIONGENERATOR_H

#include "TreeBuilder.h"
#include "OperatorAdaptor.h"
#include "CrossCorrelationCalculator.h"
#include "OperatorTree.h"
#include "Timer.h"

class CrossCorrelationGenerator : public TreeBuilder<2> {
public:
    CrossCorrelationGenerator(const MultiResolutionAnalysis<2> &mra,
                              double prec)
            : TreeBuilder<2>(mra, 30),
              build_prec(prec) { }
    virtual ~CrossCorrelationGenerator() { }

    OperatorTree *operator()(FunctionTree<1> &inp) {
        OperatorTree *out = new OperatorTree(this->MRA, this->build_prec);
        (*this)(*out, inp);
        return out;
    }

    void operator()(OperatorTree &out, FunctionTree<1> &inp) {
        this->adaptor = new OperatorAdaptor(this->build_prec, this->MRA.getMaxScale());
        this->calculator = new CrossCorrelationCalculator(inp, this->build_prec);
        this->build(out);
        this->clearCalculator();
        this->clearAdaptor();

        Timer trans_t;
        trans_t.restart();
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        out.setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
protected:
    const double build_prec;
};

#endif // CROSSCORRELATIONGENERATOR_H
