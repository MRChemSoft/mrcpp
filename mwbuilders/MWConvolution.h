#ifndef MWCONVOLUTION_H
#define MWCONVOLUTION_H

#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "OperApplicationCalculator.h"
#include "ConvolutionOperator.h"

template<int D>
class MWConvolution {
public:
    MWConvolution(ConvolutionOperator<D> &op, double pr = -1.0, int ms = MaxScale)
        : prec(pr), maxScale(ms), oper(&op) { }
    virtual ~MWConvolution() { }

    double getPrecision() const { return this->prec; }
    int getMaxScale() const { return this->maxScale; }

    void setPrecision(double pr) { this->prec = pr; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int max_iter = -1) const {
        Timer pre_t;
        this->oper->calcBandWidths(this->prec);
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        OperApplicationCalculator<D> calculator(0, this->prec, *this->oper, inp);
        pre_t.stop();

        TreeBuilder<D> builder;
        builder.build(out, calculator, adaptor, max_iter);

        Timer post_t;
        this->oper->clearBandWidths();
        out.mwTransform(TopDown, false); // add coarse scale contributions
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        inp.deleteGenerated();
        post_t.stop();

        println(10, "Time pre operator   " << pre_t);
        println(10, "Time post operator  " << post_t);
        println(10, std::endl);
    }
protected:
    double prec;
    int maxScale;
    ConvolutionOperator<D> *oper;
};

#endif // MWCONVOLUTION_H
