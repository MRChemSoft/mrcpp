#ifndef MWFDOPERATOR_H
#define MWFDOPERATOR_H

#include "TreeBuilder.h"

template<int D>
class MWFDOperator {
public:
    MWFDOperator(const MultiResolutionAnalysis<D> &mra,
                 int k, int n = -1, double pr = -1.0);
    virtual ~MWFDOperator() { }

    void setPrecision(double pr) { this->prec = pr; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int maxIter = -1,
                    int dir = -1) const;

protected:
    double prec;
    int maxScale;
    int diff_order;
    int approx_order;
};

#endif // MWFDOPERATOR_H
