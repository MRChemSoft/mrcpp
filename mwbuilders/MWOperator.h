#ifndef MWOPERATOR_H
#define MWOPERATOR_H

#include "TreeBuilder.h"
#include "OperatorTreeVector.h"

template<int D>
class MWOperator : public TreeBuilder<D> {
public:
    MWOperator(const MultiResolutionAnalysis<D> &mra, double pr)
        : TreeBuilder<D>(mra),
          apply_prec(pr),
          apply_dir(-1) {
    }

    virtual ~MWOperator() {
    }

    void setPrecision(double pr) { this->apply_prec = pr; }
    void multPrecision(double fac) { this->apply_prec *= fac; }
    void setApplyDir(int dir) { this->apply_dir = dir; }

    int getNTerms() const { return this->oper.size(); }
    const OperatorTree &getComponent(int i) const { return this->oper.getComponent(i); }

    FunctionTree<D> *operator()(FunctionTree<D> &inp);
    void operator()(FunctionTree<D> &out, FunctionTree<D> &inp, int maxIter = -1);

protected:
    int apply_dir;
    double apply_prec;
    OperatorTreeVector oper;

    void clearOperator();
};

#endif // MWOPERATOR_H
