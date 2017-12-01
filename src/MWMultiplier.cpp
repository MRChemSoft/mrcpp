#include "MWMultiplier.h"
#include "TreeBuilder.h"
#include "FunctionTree.h"
#include "WaveletAdaptor.h"
#include "MultiplicationCalculator.h"
#include "Printer.h"
#include "Timer.h"

template<int D>
void MWMultiplier<D>::operator()(FunctionTree<D> &out, double c,
                                 FunctionTree<D> &tree_a,
                                 FunctionTree<D> &tree_b,
                                 int maxIter) const {
    FunctionTreeVector<D> tree_vec;
    tree_vec.push_back(c, &tree_a);
    tree_vec.push_back(1.0, &tree_b);
    (*this)(out, tree_vec, maxIter);
}

template<int D>
void MWMultiplier<D>::operator()(FunctionTree<D> &out,
                                 FunctionTreeVector<D> &inp,
                                 int maxIter) const {
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
    MultiplicationCalculator<D> calculator(inp);

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

template class MWMultiplier<1>;
template class MWMultiplier<2>;
template class MWMultiplier<3>;
