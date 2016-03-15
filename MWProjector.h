#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

#include "TreeBuilder.h"
#include "ProjectionCalculator.h"
#include "FunctionTree.h"
#include "Timer.h"

template<int D>
class MWProjector : public TreeBuilder<D> {
public:
    MWProjector(const MultiResolutionAnalysis<D> &mra)
            : TreeBuilder<D>(mra, -1) {
        this->adaptor = new TreeAdaptor<D>();
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
};

#endif // MWPROJECTOR_H
