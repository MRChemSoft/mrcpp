#ifndef DERIVATIVEGENERATOR_H
#define DERIVATIVEGENERATOR_H

#include "TreeBuilder.h"
#include "BandWidthAdaptor.h"
#include "DerivativeCalculator.h"
#include "OperatorTree.h"
#include "Timer.h"

class DerivativeGenerator : public TreeBuilder<2> {
public:
    DerivativeGenerator(const ScalingBasis &bas, int max_scale = MaxScale)
        : TreeBuilder<2>(-1.0, max_scale), basis(bas) { }
    virtual ~DerivativeGenerator() { }

    void operator()(OperatorTree &out, double a, double b, int maxIter = -1) {
        int bw = 0;
        if (fabs(a) > MachineZero) bw = 1;
        if (fabs(b) > MachineZero) bw = 1;

        DerivativeCalculator calculator(this->basis, a, b);
        BandWidthAdaptor adaptor(bw, this->maxScale);
        this->build(out, calculator, adaptor, maxIter);

        Timer trans_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        out.setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
protected:
    ScalingBasis basis;
};

#endif // DERIVATIVEGENERATOR_H
