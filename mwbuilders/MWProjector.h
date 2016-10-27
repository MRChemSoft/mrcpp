#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

#include "TreeBuilder.h"
#include "ProjectionCalculator.h"
#include "AnalyticFunction.h"
#include "WaveletAdaptor.h"
#include "Timer.h"

template<int D>
class MWProjector : public TreeBuilder<D> {
public:
    MWProjector(double pr = -1.0, int max_scale = MaxScale)
        : TreeBuilder<D>(pr, max_scale) { }
    virtual ~MWProjector() { }

    void operator()(FunctionTree<D> &out,
                    std::function<double (const double *r)> func,
                    int maxIter = -1) {
        AnalyticFunction<D> inp(func);
        (*this)(out, inp, maxIter);
    }
    void operator()(FunctionTree<D> &out,
                    RepresentableFunction<D> &inp,
                    int maxIter = -1) {
        ProjectionCalculator<D> calculator(inp);
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        this->build(out, calculator, adaptor, maxIter);

        Timer trans_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
};

#endif // MWPROJECTOR_H
