#ifndef MWFDOPERATOR_H
#define MWFDOPERATOR_H

#include "TreeBuilder.h"

template<int D>
class MWFDOperator : public TreeBuilder<D> {
public:
    MWFDOperator(const MultiResolutionAnalysis<D> &mra,
                 int k,
                 int n = -1,
                 double prec = -1.0);
    virtual ~MWFDOperator() { }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int maxIter = -1) const;

protected:
    int diff_order;
    int approx_order;
};

#endif // MWFDOPERATOR_H
