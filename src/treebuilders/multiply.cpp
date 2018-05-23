#include "multiply.h"
#include "add.h"
#include "grid.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "MultiplicationCalculator.h"
#include "trees/HilbertIterator.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionNode.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

template<int D>
void multiply(double prec,
              FunctionTree<D> &out,
              double c,
              FunctionTree<D> &inp_a,
              FunctionTree<D> &inp_b,
              int maxIter) {
    FunctionTreeVector<D> tmp_vec;
    tmp_vec.push_back(std::make_tuple(c, &inp_a));
    tmp_vec.push_back(std::make_tuple(1.0, &inp_b));
    multiply(prec, out, tmp_vec, maxIter);
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
        FunctionTree<D> &tree = get_func(inp, i);
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

template<int D>
void dot(double prec, FunctionTree<D> &out, FunctionTreeVector<D> &inp_a, FunctionTreeVector<D> &inp_b, int maxIter) {
    if (inp_a.size() != inp_b.size()) MSG_FATAL("Input length mismatch");

    FunctionTreeVector<D> tmp_vec;
    for (int d = 0; d < inp_a.size(); d++) {
        double coef_a = get_coef(inp_a, d);
        double coef_b = get_coef(inp_b, d);
        FunctionTree<D> &tree_a = get_func(inp_a, d);
        FunctionTree<D> &tree_b = get_func(inp_b, d);
        if (out.getMRA() != tree_a.getMRA()) MSG_FATAL("Trees not compatible");
        if (out.getMRA() != tree_b.getMRA()) MSG_FATAL("Trees not compatible");
        FunctionTree<D> *out_d = new FunctionTree<D>(out.getMRA());
        copy_grid(*out_d, out);
        multiply(prec, *out_d, 1.0, tree_a, tree_b, maxIter);
        tmp_vec.push_back(std::make_tuple(coef_a*coef_b, out_d));
    }
    copy_grid(out, tmp_vec);
    add(-1.0, out, tmp_vec, 0);
    clear(tmp_vec, true);
}

template<int D>
double dot(FunctionTree<D> &bra, FunctionTree<D> &ket) {
    if (bra.getMRA() != ket.getMRA()){
        MSG_FATAL("Trees not compatible");
    }
    MWNodeVector nodeTable;
    HilbertIterator<D> it(&bra);
    it.setReturnGenNodes(false);
    while(it.next()) {
        MWNode<D> &node = it.getNode();
        nodeTable.push_back(&node);
    }
    int nNodes = nodeTable.size();
    double result = 0.0;
    double locResult = 0.0;
//OMP is disabled in order to get EXACT results (to the very last digit), the
//order of summation makes the result different beyond the 14th digit or so.
//OMP does improve the performace, but its not worth it for the time being.
//#pragma omp parallel firstprivate(n_nodes, locResult)
//		shared(nodeTable,rhs,result)
//    {
//#pragma omp for schedule(guided)
    for (int n = 0; n < nNodes; n++) {
        const FunctionNode<D> &braNode = static_cast<const FunctionNode<D> &>(*nodeTable[n]);
        const MWNode<D> *mwNode = ket.findNode(braNode.getNodeIndex());
        if (mwNode == 0) continue;

        const FunctionNode<D> &ketNode = static_cast<const FunctionNode<D> &>(*mwNode);
        if (braNode.isRootNode()) {
            locResult += dotScaling(braNode, ketNode);
        }
        locResult += dotWavelet(braNode, ketNode);
    }
//#pragma omp critical
    result += locResult;
//    }
    return result;
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
template void dot(double prec, FunctionTree<1> &out, FunctionTreeVector<1> &inp_a, FunctionTreeVector<1> &inp_b, int maxIter);
template void dot(double prec, FunctionTree<2> &out, FunctionTreeVector<2> &inp_a, FunctionTreeVector<2> &inp_b, int maxIter);
template void dot(double prec, FunctionTree<3> &out, FunctionTreeVector<3> &inp_a, FunctionTreeVector<3> &inp_b, int maxIter);
template double dot(FunctionTree<1> &bra, FunctionTree<1> &ket);
template double dot(FunctionTree<2> &bra, FunctionTree<2> &ket);
template double dot(FunctionTree<3> &bra, FunctionTree<3> &ket);

} //namespace mrcpp
