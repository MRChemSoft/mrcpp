#ifndef MWMULTIPLIER_H
#define MWMULTIPLIER_H

#include "TreeBuilder.h"
#include "MultiplicationCalculator.h"
#include "WaveletAdaptor.h"
#include "Timer.h"

template<int D>
class MWMultiplier : public TreeBuilder<D> {
public:
    MWMultiplier(double pr = -1.0) : prec(pr) { }
    virtual ~MWMultiplier() { }

    double getPrecision() const { return this->prec; }
    void setPrecision(double pr) { this->prec = pr; }
    void multPrecision(double fac) { this->prec *= fac; }

    void operator()(FunctionTree<3> &out, double c,
                    FunctionTree<D> &tree_a,
                    FunctionTree<D> &tree_b,
                    int maxIter = -1) {
        FunctionTreeVector<D> tree_vec;
        tree_vec.push_back(c, &tree_a);
        tree_vec.push_back(1.0, &tree_b);
        return (*this)(out, tree_vec, maxIter);
    }
    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int maxIter = -1) {
        this->adaptor = new WaveletAdaptor<D>(this->prec, MaxScale);
        this->calculator = new MultiplicationCalculator<D>(inp);
        this->build(out, maxIter);
        this->clearCalculator();
        this->clearAdaptor();

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

protected:
    double prec;
};

#endif // MWMULTIPLIER_H
