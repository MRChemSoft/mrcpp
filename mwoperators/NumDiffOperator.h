#ifndef NUMDIFFOPERATOR_H
#define NUMDIFFOPERATOR_H

#include "MWOperator.h"
#include "TreeBuilder.h"
#include "NumDiffCalculator.h"
#include "BandWidthAdaptor.h"

template<int D>
class NumDiffOperator : public MWOperator {
public:
    NumDiffOperator(const MultiResolutionAnalysis<D> &mra)
            : MWOperator(mra.getOperatorMRA()) {
        initializeOperator();
    }
    virtual ~NumDiffOperator() { this->clearOperator(); }
    bool applyCompressed() const { return false; }

protected:
    void initializeOperator() {
        int bw = 1; //Operator bandwidth
        int max_scale = this->oper_mra.getMaxScale();
        const ScalingBasis &basis = this->oper_mra.getScalingBasis();
        NumDiffCalculator calculator(basis);
        BandWidthAdaptor adaptor(bw, max_scale);
        TreeBuilder<2> builder;

        OperatorTree *o_tree = new OperatorTree(this->oper_mra, MachineZero, MaxAllocOperNodes);
        builder.build(*o_tree, calculator, adaptor, -1);

        Timer trans_t;
        //o_tree->mwTransform(BottomUp);
        o_tree->calcSquareNorm();
        o_tree->setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);

        this->oper_exp.push_back(o_tree);
    }
};

#endif // NUMDIFFOPERATOR_H
