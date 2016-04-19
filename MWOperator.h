#ifndef MWOPERATOR_H
#define MWOPERATOR_H

#include "TreeBuilder.h"
#include "FunctionTree.h"
#include "OperApplicationCalculator.h"
#include "OperatorTreeVector.h"
#include "Timer.h"

template<int D>
class MWOperator : public TreeBuilder<D> {
public:
    MWOperator(const MultiResolutionAnalysis<D> &mra,
               double prec = -1.0, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->adaptor = new WaveletAdaptor<D>(prec);
    }
    MWOperator(const MultiResolutionAnalysis<D> &mra,
               const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
       this->adaptor = a.copy();
    }
    virtual ~MWOperator() {
        this->clearAdaptor();
    }

    FunctionTree<D> *operator()(FunctionTree<D> &inp) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        initializeGrid(*out, inp);
        (*this)(*out, inp);
        return out;
    }

    void operator()(FunctionTree<D> &out, FunctionTree<D> &inp) {
        this->calculator = new OperApplicationCalculator<D>(this->oper, inp);
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
    OperatorTreeVector oper;

    /** Build grid based on analytic input function */
    void initializeGrid(FunctionTree<D> &out, FunctionTree<D> &inp) {
        Timer init_t;
        init_t.restart();
        GridGenerator<D> G(this->MRA);
        G(out, inp);
        init_t.stop();
        println(10, "Time initializing   " << init_t);
    }
};

#endif // MWOPERATOR_H
