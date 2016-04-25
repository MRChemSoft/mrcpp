#ifndef MWOPERATOR_H
#define MWOPERATOR_H

#include "TreeBuilder.h"
#include "OperatorTreeVector.h"

template<int D>
class MWOperator : public TreeBuilder<D> {
public:
    MWOperator(const MultiResolutionAnalysis<D> &mra, double prec, int iter);
    virtual ~MWOperator();

    void setApplyPrecision(double prec) { this->apply_prec = prec; }
    void multApplyPrecision(double fac) { this->apply_prec *= fac; }

    FunctionTree<D> *operator()(FunctionTree<D> &inp);
    void operator()(FunctionTree<D> &out, FunctionTree<D> &inp);

protected:
    double apply_prec;
    OperatorTreeVector oper;

    void clearOperator();
    void initializeGrid(FunctionTree<D> &out, FunctionTree<D> &inp);
    MultiResolutionAnalysis<2> *getOperatorMRA();
};

#endif // MWOPERATOR_H
