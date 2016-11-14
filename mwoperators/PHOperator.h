#ifndef PHOPERATOR_H
#define PHOPERATOR_H

#include "DerivativeOperator.h"
#include "TreeBuilder.h"
#include "PHCalculator.h"
#include "BandWidthAdaptor.h"

template<int D>
class PHOperator : public DerivativeOperator {
public:
    PHOperator(const MultiResolutionAnalysis<D> &mra, int order)
            : MWOperator(mra.getOperatorMRA()) {
        initializeOperator(order);
    }
    virtual ~PHOperator() { }

protected:
    void initializeOperator(int order) {
        int bw = 1; //Operator bandwidth
        int max_scale = this->oper_mra.getMaxScale();
        const ScalingBasis &basis = this->oper_mra.getScalingBasis();
        PHCalculator calculator(basis, order);
        BandWidthAdaptor adaptor(bw, max_scale);
        TreeBuilder<2> builder;
        builder.build(this->oper_tree, calculator, adaptor, -1);

        Timer trans_t;
        //o_tree->mwTransform(BottomUp);
        this->oper_tree.calcSquareNorm();
        this->oper_tree.setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
};

#endif // PHOPERATOR_H
