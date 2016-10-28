#ifndef OPERATORAPPLIER_H
#define OPERATORAPPLIER_H

#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "OperApplicationCalculator.h"
#include "MWOperator.h"

template<int D>
class OperatorApplier : public TreeBuilder<D> {
public:
    OperatorApplier(double pr = -1.0, int max_scale = MaxScale)
        : TreeBuilder<D>(pr, max_scale) { }
    virtual ~OperatorApplier() { }

    void operator()(FunctionTree<D> &out,
                    MWOperator &oper,
                    FunctionTree<D> &inp,
                    int maxIter = -1,
                    int dir = -1) const {
        Timer pre_t;
        oper.calcBandWidths(this->prec);
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        OperApplicationCalculator<D> calculator(dir, this->prec, oper, inp);
        pre_t.stop();

        this->build(out, calculator, adaptor, maxIter);

        Timer post_t;
        oper.clearBandWidths();
        out.mwTransform(TopDown, false); // add coarse scale contributions
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        inp.deleteGenerated();
        post_t.stop();

        println(10, "Time pre operator   " << pre_t);
        println(10, "Time post operator  " << post_t);
        println(10, std::endl);
    }
};

#endif // OPERATORAPPLIER_H
