#include "add.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "AdditionCalculator.h"
#include "FunctionTree.h"
#include "Printer.h"
#include "Timer.h"

using namespace mrcpp;

template<int D>
void mrcpp::add(double prec, FunctionTree<D> &out,
                double a, FunctionTree<D> &tree_a,
                double b, FunctionTree<D> &tree_b,
                int maxIter) {
    FunctionTreeVector<D> tree_vec;
    tree_vec.push_back(a, &tree_a);
    tree_vec.push_back(b, &tree_b);
    mrcpp::add(prec, out, tree_vec, maxIter);
}

template<int D>
void mrcpp::add(double prec,
                FunctionTree<D> &out,
                FunctionTreeVector<D> &inp,
                int maxIter) {
    int maxScale = out.getMRA().getMaxScale();
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(prec, maxScale);
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

    Printer::printTime(10, "Time transform", trans_t);
    Printer::printTime(10, "Time cleaning", clean_t);
    Printer::printSeparator(10, ' ');
}

template void mrcpp::add(double prec, FunctionTree<1> &out, double a, FunctionTree<1> &tree_a, double b, FunctionTree<1> &tree_b, int maxIter);
template void mrcpp::add(double prec, FunctionTree<2> &out, double a, FunctionTree<2> &tree_a, double b, FunctionTree<2> &tree_b, int maxIter);
template void mrcpp::add(double prec, FunctionTree<3> &out, double a, FunctionTree<3> &tree_a, double b, FunctionTree<3> &tree_b, int maxIter);
template void mrcpp::add(double prec, FunctionTree<1> &out, FunctionTreeVector<1> &inp, int maxIter);
template void mrcpp::add(double prec, FunctionTree<2> &out, FunctionTreeVector<2> &inp, int maxIter);
template void mrcpp::add(double prec, FunctionTree<3> &out, FunctionTreeVector<3> &inp, int maxIter);
