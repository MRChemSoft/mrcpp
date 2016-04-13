#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

#include "TreeBuilder.h"
#include "ProjectionCalculator.h"
#include "FunctionTree.h"
#include "Timer.h"

template<int D>
class MWProjector : public TreeBuilder<D> {
public:
    MWProjector(const MultiResolutionAnalysis<D> &mra,
                double prec = -1.0, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->adaptor = new WaveletAdaptor<D>(prec);
    }
    MWProjector(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
       this->adaptor = a.copy();
    }
    virtual ~MWProjector() {
        this->clearAdaptor();
    }

    FunctionTree<D> *operator()(RepresentableFunction<D> &inp) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        initializeGrid(*out, inp);
        (*this)(*out, inp);
        return out;
    }

    void operator()(FunctionTree<D> &out, RepresentableFunction<D> &inp) {
        this->calculator = new ProjectionCalculator<D>(inp);
        this->build(out);
        this->clearCalculator();

        Timer trans_t;
        trans_t.restart();
        out.mwTransform(BottomUp);
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);
    }
protected:
    /** Build grid based on analytic input function */
    void initializeGrid(FunctionTree<D> &out, RepresentableFunction<D> &inp) {
        Timer init_t;
        init_t.restart();
        GridGenerator<D> G(this->MRA);
        G(out, inp);
        init_t.stop();
        println(10, "Time initializing   " << init_t);
    }
};

#endif // MWPROJECTOR_H
