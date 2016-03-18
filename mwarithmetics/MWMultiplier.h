#ifndef MWMULTIPLIER_H
#define MWMULTIPLIER_H

#include "TreeBuilder.h"
#include "MultiplicationCalculator.h"
#include "GridGenerator.h"

template<int D>
class MWMultiplier : public TreeBuilder<D> {
public:
    MWMultiplier(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->adaptor = new TreeAdaptor<D>();
    }
    MWMultiplier(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->adaptor = a.copy();
    }
    virtual ~MWMultiplier() {
        this->clearAdaptor();
    }

    FunctionTree<D>* operator()(MultiplicationVector<D> &inp) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        initializeGrid(*out, inp);
        (*this)(*out, inp);
        return out;
    }

    void operator()(FunctionTree<D> &out, MultiplicationVector<D> &inp) {
        Timer trans_t, clean_t;
        this->calculator = new MultiplicationCalculator<D>(inp);
        this->build(out);
        this->clearCalculator();

        trans_t.restart();
        out.mwTransform(BottomUp);
        trans_t.stop();

        clean_t.restart();
        inp.clean();
        clean_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, "Time cleaning       " << clean_t);
        println(10, std::endl);
    }
protected:
    /** Copy the grids from the input functions */
    void initializeGrid(FunctionTree<D> &out, MultiplicationVector<D> &inp) {
        Timer init_t;
        init_t.restart();
        GridGenerator<D> G(this->MRA);
        for (int i = 0; i < inp.size(); i++) {
            FunctionTree<D> &func_i = inp.getFunc(i);
            G(out, func_i);
        }
        init_t.stop();
        println(10, "Time initializing   " << init_t);
    }
};


#endif // MWMULTIPLIER_H
