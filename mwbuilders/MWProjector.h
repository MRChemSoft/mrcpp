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
    MWProjector(double pr = -1.0) : prec(pr) { }
    virtual ~MWProjector() { }

    void setPrecision(double pr) { this->prec = pr; }
    void multPrecision(double fac) { this->prec *= fac; }

    void operator()(FunctionTree<D> &out,
                    std::function<double (const double *r)> func,
                    int maxIter = -1) {
        AnalyticFunction<D> inp(func);
        (*this)(out, inp, maxIter);
    }

    void operator()(FunctionTree<D> &out,
                    RepresentableFunction<D> &inp,
                    int maxIter = -1) {
        this->adaptor = new WaveletAdaptor<D>(this->prec, MaxScale);
        this->calculator = new ProjectionCalculator<D>(inp);
        this->build(out, maxIter);
        this->clearCalculator();
        this->clearAdaptor();

        Timer trans_t;
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
protected:
    double prec;
};

#endif // MWPROJECTOR_H
