#ifndef MWADDER_H
#define MWADDER_H

#include "TreeBuilder.h"
#include "AdditionCalculator.h"
#include "WaveletAdaptor.h"
#include "Timer.h"
#include "TreeAllocator.h"

template<int D>
class MWAdder : public TreeBuilder<D> {
public:
    MWAdder(const MultiResolutionAnalysis<D> &mra, double pr = -1.0)
            : TreeBuilder<D>(mra),
              prec(pr) {
    }
    virtual ~MWAdder() {
    }

    double getPrecision() const { return this->prec; }
    void setPrecision(double pr) { this->prec = pr; }
    void multPrecision(double fac) { this->prec *= fac; }

    FunctionTree<D>* operator()(double a, FunctionTree<D> &tree_a,
                                double b, FunctionTree<D> &tree_b) {
        FunctionTreeVector<D> tree_vec;
        tree_vec.push_back(a, &tree_a);
        tree_vec.push_back(b, &tree_b);
        return (*this)(tree_vec);
    }
    FunctionTree<D>* operator()(FunctionTreeVector<D> &inp) {
      TreeAllocator<D> *alloc = new TreeAllocator<D>(this->MRA, MAXALLOCNODES);
      FunctionTree<D> *out = alloc->getTree();
       //       FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp);
        return out;
    }

    void operator()(FunctionTree<3> &out,
                    double a, FunctionTree<D> &tree_a,
                    double b, FunctionTree<D> &tree_b,
                    int maxIter = -1) {
        FunctionTreeVector<D> tree_vec;
        tree_vec.push_back(a, &tree_a);
        tree_vec.push_back(b, &tree_b);
        return (*this)(out, tree_vec, maxIter);
    }
    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int maxIter = -1) {
        this->adaptor = new WaveletAdaptor<D>(this->prec, this->MRA.getMaxScale());
        this->calculator = new AdditionCalculator<D>(inp);
        this->build(out, maxIter);
        this->clearCalculator();
        this->clearAdaptor();

        Timer trans_t;
        trans_t.restart();
        out.mwTransform(BottomUp);
        out.calcSquareNorm();
        trans_t.stop();

        Timer clean_t;
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
    double prec;
};

#endif // MWADDER_H
