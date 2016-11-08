#ifndef HARRISONDERIVATIVE_H
#define HARRISONDERIVATIVE_H

#include "TreeBuilder.h"
#include "HarrisonDerivativeCalculator.h"
#include "WaveletAdaptor.h"
#include "FunctionTree.h"
#include "Timer.h"

template<int D>
class HarrisonDerivative : public TreeBuilder<D> {
public:
    HarrisonDerivative(double prec = -1.0) 
        : TreeBuilder<D>(prec, MaxScale) { }
    virtual ~HarrisonDerivative() { }

    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int max_iter = -1,
                    int dir = -1) const {
        Timer pre_t;
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        HarrisonDerivativeCalculator<D> calculator(dir, inp);
        pre_t.stop();

        this->build(out, calculator, adaptor, max_iter);

        Timer post_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        inp.deleteGenerated();
        post_t.stop();

        println(10, "Time pre operator   " << pre_t);
        println(10, "Time post operator  " << post_t);
        println(10, std::endl);
    }
};

#endif // HARRISONDERIVATIVE_H
