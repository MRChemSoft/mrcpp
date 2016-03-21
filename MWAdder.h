#ifndef MWADDER_H
#define MWADDER_H

#include "TreeBuilder.h"
#include "AdditionCalculator.h"
#include "GridGenerator.h"

template<int D>
class MWAdder : public TreeBuilder<D> {
public:
    MWAdder(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->adaptor = new TreeAdaptor<D>();
    }
    MWAdder(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->adaptor = a.copy();
    }
    virtual ~MWAdder() {
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
        this->calculator = new AdditionCalculator<D>(inp);
        this->build(out);
        this->clearCalculator();

        trans_t.restart();
        out.mwTransform(BottomUp);
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

#endif // MWADDER_H
