#ifndef CROSSCORRELATIONGENERATOR_H
#define CROSSCORRELATIONGENERATOR_H

#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "CrossCorrelationCalculator.h"
#include "OperatorTree.h"

class CrossCorrelationGenerator : public TreeBuilder<2> {
public:
    CrossCorrelationGenerator(const MultiResolutionAnalysis<2> &mra,
                              double prec)
            : TreeBuilder<2>(mra, 30),
              build_prec(prec) { }
    virtual ~CrossCorrelationGenerator() { }

    OperatorTree *operator()(FunctionTree<1> &kernel) {
        this->adaptor = new WaveletAdaptor<2>(this->build_prec);
        this->calculator = new CrossCorrelationCalculator(kernel, this->build_prec);
        OperatorTree *out = new OperatorTree(this->MRA);
        this->build(*out);
        out->mwTransform(BottomUp);
        out->setupOperNodeCache();
        this->clearCalculator();
        this->clearAdaptor();
        return out;
    }
protected:
    const double build_prec;
};

#endif // CROSSCORRELATIONGENERATOR_H
