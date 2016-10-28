#ifndef MWOPERATOR_H
#define MWOPERATOR_H

#include <vector>

#include "MultiResolutionAnalysis.h"
#include "OperatorTree.h"

class MWOperator {
public:
    MWOperator(MultiResolutionAnalysis<2> mra) : oper_mra(mra) { }
    virtual ~MWOperator() { }

    int size() const { return this->oper_exp.size(); }
    void push_back(OperatorTree *oper) { this->oper_exp.push_back(oper); } 
    void clear(bool dealloc = false) {
        if (dealloc) {
            for (int i = 0; i < oper_exp.size(); i++) {
                if (this->oper_exp[i] != 0) delete this->oper_exp[i];
            }
        }
        this->oper_exp.clear();
    }

    int getMaxBandWidth(int depth = -1) const;
    const Eigen::VectorXi &getMaxBandWidths() const { return this->bandMax; }

    void calcBandWidths(double prec);
    void clearBandWidths();

    OperatorTree &getComponent(int i);
    const OperatorTree &getComponent(int i) const;

    OperatorTree *operator[](int i) { return this->oper_exp[i]; }
    const OperatorTree *operator[](int i) const { return this->oper_exp[i]; }

protected:
    MultiResolutionAnalysis<2> oper_mra;
    std::vector<OperatorTree *> oper_exp;
    Eigen::VectorXi bandMax;

    void clearOperator();
};

#endif // MWOPERATOR_H
