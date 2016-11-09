#ifndef DERIVATIVEOPERATOR_H
#define DERIVATIVEOPERATOR_H

#include "MWOperator.h"
#include "TreeBuilder.h"
#include "DerivativeCalculator.h"
#include "BandWidthAdaptor.h"

template<int D>
class DerivativeOperator : public MWOperator {
public:
    DerivativeOperator(const MultiResolutionAnalysis<D> &mra,
                       double a = 0.5, double b = 0.5)
            : MWOperator(mra.getOperatorMRA()) {
        initializeOperator(a, b);
    }
    virtual ~DerivativeOperator() { this->clearOperator(); }

protected:
    void initializeOperator(double a, double b) {
        int bw = 0; //Operator bandwidth
        if (fabs(a) > MachineZero) bw = 1;
        if (fabs(b) > MachineZero) bw = 1;

        int max_scale = this->oper_mra.getMaxScale();
        const ScalingBasis &basis = this->oper_mra.getScalingBasis();
        DerivativeCalculator calculator(basis, a, b);
        BandWidthAdaptor adaptor(bw, max_scale);
        TreeBuilder<2> builder;

        OperatorTree *o_tree = new OperatorTree(this->oper_mra, MachineZero, MaxAllocOperNodes);
        builder.build(*o_tree, calculator, adaptor, -1);

        Timer trans_t;
        o_tree->mwTransform(BottomUp);
        o_tree->calcSquareNorm();
        o_tree->setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);

        this->oper_exp.push_back(o_tree);
    }
};

#endif // DERIVATIVEOPERATOR_H
