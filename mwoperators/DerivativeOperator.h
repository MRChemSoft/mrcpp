#ifndef DERIVATIVEOPERATOR_H
#define DERIVATIVEOPERATOR_H

#include "MultiResolutionAnalysis.h"
#include "OperatorTree.h"

class DerivativeOperator {
public:
    DerivativeOperator(MultiResolutionAnalysis<2> mra) 
        : oper_tree(mra, MachineZero, MaxAllocOperNodes) { }
    virtual ~DerivativeOperator() { }

    OperatorTree &getOperatorTree() { return this->oper_tree; }
    const OperatorTree &getOperatorTree() const { return this->oper_tree; }

protected:
    OperatorTree oper_tree;
};

#endif // DERIVATIVEOPERATOR_H
