#include "MWAdder.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "AdditionCalculator.h"
#include "FunctionTree.h"
#include "Printer.h"
#include "Timer.h"

template<int D>
void MWAdder<D>::operator()(FunctionTree<D> &out,
                            double a, FunctionTree<D> &tree_a,
                            double b, FunctionTree<D> &tree_b,
                            int maxIter) const {
    FunctionTreeVector<D> tree_vec;
    tree_vec.push_back(a, &tree_a);
    tree_vec.push_back(b, &tree_b);
    return (*this)(out, tree_vec, maxIter);
}

template<int D>
void MWAdder<D>::operator()(FunctionTree<D> &out,
                            FunctionTreeVector<D> &inp,
                            int maxIter) const {
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
    AdditionCalculator<D> calculator(inp);

    builder.build(out, calculator, adaptor, maxIter);

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

template class MWAdder<1>;
template class MWAdder<2>;
template class MWAdder<3>;
