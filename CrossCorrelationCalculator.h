#ifndef CROSSCORRELATIONCALCULATOR_H
#define CROSSCORRELATIONCALCULATOR_H

#include "TreeCalculator.h"
#include "CrossCorrelationCache.h"

class CrossCorrelationCalculator : public TreeCalculator<2> {
public:
    CrossCorrelationCalculator(FunctionTree<1> &k, double p)
            : prec(p),
              kernel(&k) { }
    virtual ~CrossCorrelationCalculator() { }

protected:
    const double prec;
    FunctionTree<1> *kernel;

    virtual void calcNode(MWNode<2> &node);

    template<int T>
    void applyCcc(MWNode<2> &node, CrossCorrelationCache<T> &ccc);
};


#endif // CROSSCORRELATIONCALCULATOR_H
