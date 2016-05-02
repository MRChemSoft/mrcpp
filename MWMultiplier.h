#ifndef MWMULTIPLIER_H
#define MWMULTIPLIER_H

#include "TreeBuilder.h"
#include "MultiplicationCalculator.h"
#include "GridGenerator.h"

template<int D>
class MWMultiplier : public TreeBuilder<D> {
public:
    MWMultiplier(const MultiResolutionAnalysis<D> &mra,
                 double prec = -1.0,int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->adaptor = new WaveletAdaptor<D>(prec, mra.getMaxScale());
    }
    MWMultiplier(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->adaptor = a.copy();
    }
    virtual ~MWMultiplier() {
        this->clearAdaptor();
    }

    FunctionTree<D>* operator()(FunctionTreeVector<D> &inp) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        initializeGrid(*out, inp);
        (*this)(*out, inp);
        return out;
    }

    void operator()(FunctionTree<D> &out, FunctionTreeVector<D> &inp) {
        Timer trans_t, clean_t;
        this->calculator = new MultiplicationCalculator<D>(inp);
        this->build(out);
        this->clearCalculator();

        trans_t.restart();
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        trans_t.stop();

        clean_t.restart();
        for (int i = 0; i < inp.size(); i++) {
            FunctionTree<D> &tree = inp.getFunc(i);
            tree.deleteGenerated();
        }
        clean_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, "Time cleaning       " << clean_t);
        println(10, std::endl);
    }
protected:
    /** Copy the grids from the input functions */
    void initializeGrid(FunctionTree<D> &out, FunctionTreeVector<D> &inp) {
        Timer init_t;
        init_t.restart();
        GridGenerator<D> G(this->MRA);
        G(out, inp);
        init_t.stop();
        println(10, "Time initializing   " << init_t);
    }
};


#endif // MWMULTIPLIER_H
