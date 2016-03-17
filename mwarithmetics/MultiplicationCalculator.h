#ifndef MULTIPLICATIONCALCULATOR_H
#define MULTIPLICATIONCALCULATOR_H

#include "TreeCalculator.h"
#include "MultiplicationVector.h"

template<int D> class MWMultiplier;

template<int D>
class MultiplicationCalculator : public TreeCalculator<D> {
public:
    friend class MWMultiplier<D>;
protected:
    MultiplicationCalculator(MultiplicationVector<D> &inp) : prod_vec(&inp) { }
    virtual ~MultiplicationCalculator() { }

    virtual void calcNode(MWNode<D> &node_o) const {
        NOT_IMPLEMENTED_ABORT;
    }

private:
    MultiplicationVector<D> *prod_vec;
};

#endif // MULTIPLICATIONCALCULATOR_H
