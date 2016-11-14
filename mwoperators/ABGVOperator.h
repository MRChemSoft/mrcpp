#ifndef ABGVOPERATOR_H
#define ABGVOPERATOR_H

#include "DerivativeOperator.h"
#include "TreeBuilder.h"
#include "ABGVCalculator.h"
#include "BandWidthAdaptor.h"

template<int D>
class ABGVOperator : public DerivativeOperator {
public:
    ABGVOperator(const MultiResolutionAnalysis<D> &mra, double a, double b)
            : DerivativeOperator(mra.getOperatorMRA()) {
        initializeOperator(a, b);
    }
    virtual ~ABGVOperator() { }

protected:
    void initializeOperator(double a, double b) {
        int bw = 0; //Operator bandwidth
        if (fabs(a) > MachineZero) bw = 1;
        if (fabs(b) > MachineZero) bw = 1;

        int max_scale = this->oper_mra.getMaxScale();
        const ScalingBasis &basis = this->oper_mra.getScalingBasis();
        ABGVCalculator calculator(basis, a, b);
        BandWidthAdaptor adaptor(bw, max_scale);
        TreeBuilder<2> builder;
        builder.build(this->oper_tree, calculator, adaptor, -1);

        Timer trans_t;
        //this->oper_tree.mwTransform(BottomUp);
        this->oper_tree.calcSquareNorm();
        this->oper_tree.setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
};

#endif // ABGVOPERATOR_H
