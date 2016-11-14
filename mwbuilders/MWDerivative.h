#ifndef MWDERIVATIVE_H
#define MWDERIVATIVE_H

#include "TreeBuilder.h"
#include "CopyAdaptor.h"
#include "OperApplicationCalculator.h"
#include "DerivativeOperator.h"

template<int D>
class MWDerivative {
public:
    MWDerivative(DerivativeOperator<D> &op, int ms = MaxScale)
        : maxScale(ms), oper(&op) { }
    virtual ~MWDerivative() { this->oper = 0; }

    int getMaxScale() const { return this->maxScale; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int dir = -1) const {
        TreeBuilder<D> builder;

        Timer pre_t;
        int bw[D];
        getBandWidth(dir, bw);
        CopyAdaptor<D> pre_adaptor(inp, this->maxScale, bw);
        DefaultCalculator<D> pre_calculator;
        builder.build(out, pre_calculator, pre_adaptor, -1);
        pre_t.stop();

        TreeAdaptor<D> apply_adaptor(this->maxScale);
        OperApplicationCalculator<D> apply_calculator(dir, this->prec, *this->oper, inp);
        builder.build(out, apply_calculator, apply_adaptor, 0);

        Timer post_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        inp.deleteGenerated();
        post_t.stop();

        println(10, "Time pre operator   " << pre_t);
        println(10, "Time post operator  " << post_t);
        println(10, std::endl);
    }
protected:
    int maxScale;
    DerivativeOperator<D> *oper;

    void getBandWidth(int dir, int bw[D]) const {
        if (dir >= D) MSG_FATAL("Invalid apply dir");
        for (int d = 0; d < D; d++) bw[d] = 0;
        bw[dir] = 1;
    }
};

#endif // MWDERIVATIVE_H
