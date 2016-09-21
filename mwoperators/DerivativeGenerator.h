#ifndef DERIVATIVEGENERATOR_H
#define DERIVATIVEGENERATOR_H

#include "TreeBuilder.h"
#include "BandWidthAdaptor.h"
#include "DerivativeCalculator.h"
#include "OperatorTree.h"
#include "Timer.h"

class DerivativeGenerator : public TreeBuilder<2> {
public:
    DerivativeGenerator(const MultiResolutionAnalysis<2> &mra)
        : TreeBuilder<2>(mra) { }
    virtual ~DerivativeGenerator() { }

    OperatorTree *operator()(double a, double b) {
        OperatorTree *out = new OperatorTree(this->MRA, MachineZero);
        (*this)(*out, a, b, -1);
        return out;
    }

    void operator()(OperatorTree &out, double a, double b, int maxIter = -1) {
        int bw = 0;
        if (fabs(a) > MachineZero) bw = 1;
        if (fabs(b) > MachineZero) bw = 1;

        this->adaptor = new BandWidthAdaptor(bw, this->MRA.getMaxScale());
        this->calculator = new DerivativeCalculator(this->MRA.getScalingBasis(), a, b);
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
};

#endif // DERIVATIVEGENERATOR_H
