#include "multiply.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "MultiplicationCalculator.h"
#include "trees/FunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

template<int D>
void multiply(double prec,
              FunctionTree<D> &out,
              double c,
              FunctionTree<D> &tree_a,
              FunctionTree<D> &tree_b,
              int maxIter) {
    FunctionTreeVector<D> tree_vec;
    tree_vec.push_back(std::make_tuple(c, &tree_a));
    tree_vec.push_back(std::make_tuple(1.0, &tree_b));
    mrcpp::multiply(prec, out, tree_vec, maxIter);
}

template<int D>
void multiply(double prec,
              FunctionTree<D> &out,
              FunctionTreeVector<D> &inp,
              int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
    MultiplicationCalculator<D> calculator(inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    Timer clean_t;
    for (int i = 0; i < inp.size(); i++) {
        FunctionTree<D> &tree = getFunc(inp, i);
        tree.deleteGenerated();
    }
    clean_t.stop();

    Printer::printTime(10, "Time transform", trans_t);
    Printer::printTime(10, "Time cleaning", clean_t);
    Printer::printSeparator(10, ' ');
}

template<int D>
void power(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, double pow) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void map(double prec, FunctionTree<D> &out, FunctionTree<D> &inp, RepresentableFunction<D> &func) {
    NOT_IMPLEMENTED_ABORT;
}

template void multiply(double prec, FunctionTree<1> &out, double c, FunctionTree<1> &tree_a, FunctionTree<1> &tree_b, int maxIter);
template void multiply(double prec, FunctionTree<2> &out, double c, FunctionTree<2> &tree_a, FunctionTree<2> &tree_b, int maxIter);
template void multiply(double prec, FunctionTree<3> &out, double c, FunctionTree<3> &tree_a, FunctionTree<3> &tree_b, int maxIter);
template void multiply(double prec, FunctionTree<1> &out, FunctionTreeVector<1> &inp, int maxIter);
template void multiply(double prec, FunctionTree<2> &out, FunctionTreeVector<2> &inp, int maxIter);
template void multiply(double prec, FunctionTree<3> &out, FunctionTreeVector<3> &inp, int maxIter);
template void power(double prec, FunctionTree<1> &out, FunctionTree<1> &tree, double pow);
template void power(double prec, FunctionTree<2> &out, FunctionTree<2> &tree, double pow);
template void power(double prec, FunctionTree<3> &out, FunctionTree<3> &tree, double pow);
template void map(double prec, FunctionTree<1> &out, FunctionTree<1> &inp, RepresentableFunction<1> &func);
template void map(double prec, FunctionTree<2> &out, FunctionTree<2> &inp, RepresentableFunction<2> &func);
template void map(double prec, FunctionTree<3> &out, FunctionTree<3> &inp, RepresentableFunction<3> &func);

} //namespace mrcpp
