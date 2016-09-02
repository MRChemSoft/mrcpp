#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

#include "TreeBuilder.h"
#include "ProjectionCalculator.h"
#include "AnalyticFunction.h"
#include "WaveletAdaptor.h"
#include "SerialTree.h"
#include "Timer.h"

template<int D>
class MWProjector : public TreeBuilder<D> {
public:
    MWProjector(const MultiResolutionAnalysis<D> &mra, double pr = -1.0)
            : TreeBuilder<D>(mra),
              prec(pr) {
    }
    virtual ~MWProjector() {
    }

    void setPrecision(double pr) { this->prec = pr; }
    void multPrecision(double fac) { this->prec *= fac; }

    FunctionTree<D> *operator()(std::function<double (const double *r)> func) {
        AnalyticFunction<D> inp(func);
        SerialTree<D> *alloc = new SerialTree<D>(this->MRA, MAXALLOCNODES);
        FunctionTree<D> *out = alloc->getTree();
        //FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp);
        return out;
    }

    FunctionTree<D> *operator()(RepresentableFunction<D> &inp) {
        //FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        SerialTree<D> *alloc = new SerialTree<D>(this->MRA, MAXALLOCNODES);
        FunctionTree<D> *out = alloc->getTree();
        (*this)(*out, inp);
        return out;
    }

    void operator()(FunctionTree<D> &out,
                    RepresentableFunction<D> &inp,
                    int maxIter = -1) {
        this->adaptor = new WaveletAdaptor<D>(this->prec, this->MRA.getMaxScale());
        this->calculator = new ProjectionCalculator<D>(inp);
        this->build(out, maxIter);
        this->clearCalculator();
        this->clearAdaptor();

        Timer trans_t;
        trans_t.restart();
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
