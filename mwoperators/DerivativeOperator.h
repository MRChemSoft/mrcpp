#ifndef DERIVATIVEOPERATOR_H
#define DERIVATIVEOPERATOR_H

#include "MWOperator.h"
#include "DerivativeGenerator.h"

template<int D>
class DerivativeOperator : public MWOperator {
public:
    DerivativeOperator(const MultiResolutionAnalysis<D> &mra,
                       double a = 0.5, double b = 0.5)
            : MWOperator(mra.getOperatorMRA()),
              A(a), B(b) {
        if (this->A > MachineZero) NEEDS_TESTING;
        if (this->B > MachineZero) NEEDS_TESTING;
        initializeOperator();
    }
    virtual ~DerivativeOperator() { this->clearOperator(); }

protected:
    const double A;
    const double B;

    void initializeOperator() {
        int max_scale = this->oper_mra.getMaxScale();
        const ScalingBasis &basis = this->oper_mra.getScalingBasis();
        DerivativeGenerator DG(basis, max_scale);

        OperatorTree *oper_comp = new OperatorTree(this->oper_mra, MachineZero, MaxAllocOperNodes);
        DG(*oper_comp, this->A, this->B);
        this->oper_exp.push_back(oper_comp);
    }
};

#endif // DERIVATIVEOPERATOR_H
