#ifndef MWADDER_H
#define MWADDER_H

#include "TreeBuilder.h"
#include "AdditionCalculator.h"
#include "WaveletAdaptor.h"
#include "Timer.h"

template<int D>
class MWAdder : public TreeBuilder<D> {
public:
    MWAdder(double pr = -1.0, int max_scale = MaxScale)
        : TreeBuilder<D>(pr, max_scale) { }
    virtual ~MWAdder() { }

    void operator()(FunctionTree<D> &out,
                    double a, FunctionTree<D> &tree_a,
                    double b, FunctionTree<D> &tree_b,
                    int maxIter = -1) const {
        FunctionTreeVector<D> tree_vec;
        tree_vec.push_back(a, &tree_a);
        tree_vec.push_back(b, &tree_b);
        return (*this)(out, tree_vec, maxIter);
    }
    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int maxIter = -1) const {
        AdditionCalculator<D> calculator(inp);
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        this->build(out, calculator, adaptor, maxIter);

        Timer trans_t;
        out.mwTransform(BottomUp);
	out.calcSquareNorm();
        trans_t.stop();

        Timer clean_t;
        for (int i = 0; i < inp.size(); i++) {
            FunctionTree<D> &tree = inp.getFunc(i);
            tree.deleteGenerated();
        }
        clean_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, "Time cleaning       " << clean_t);
        println(10, std::endl);
    }
};

#endif // MWADDER_H
