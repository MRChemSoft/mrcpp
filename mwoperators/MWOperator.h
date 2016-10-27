#ifndef MWOPERATOR_H
#define MWOPERATOR_H

#include <vector>

#include "MultiResolutionAnalysis.h"
#include "OperatorTree.h"

template<int D>
class MWOperator {
public:
    MWOperator(const MultiResolutionAnalysis<D> &mra) : MRA(mra) { }
    virtual ~MWOperator() { }

    int size() const { return this->oper_exp.size(); }

    int getMaxBandWidth(int depth = -1) const;
    const Eigen::VectorXi &getMaxBandWidths() const { return this->bandMax; }

    void calcBandWidths(double prec);
    void clearBandWidths();

    OperatorTree &getComponent(int i);
    const OperatorTree &getComponent(int i) const;

    OperatorTree *operator[](int i) { return this->oper_exp[i]; }
    const OperatorTree *operator[](int i) const { return this->oper_exp[i]; }

protected:
    const MultiResolutionAnalysis<D> MRA;
    std::vector<OperatorTree *> oper_exp;
    Eigen::VectorXi bandMax;

    void clearOperator();
};

#endif // MWOPERATOR_H
