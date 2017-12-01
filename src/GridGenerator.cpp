#include "GridGenerator.h"
#include "TreeBuilder.h"
#include "AnalyticAdaptor.h"
#include "CopyAdaptor.h"
#include "DefaultCalculator.h"
#include "Printer.h"

template<int D>
void GridGenerator<D>::operator()(FunctionTree<D> &out,
                                  const RepresentableFunction<D> &inp,
                                  int maxIter) const {
    TreeBuilder<D> builder;
    AnalyticAdaptor<D> adaptor(inp, this->maxScale);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, maxIter);
    println(10, std::endl);
}

template<int D>
void GridGenerator<D>::operator()(FunctionTree<D> &out,
                                  FunctionTree<D> &inp,
                                  int max_iter) const {
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, this->maxScale, 0);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, max_iter);
    println(10, std::endl);
}

template<int D>
void GridGenerator<D>::operator()(FunctionTree<D> &out,
                                  FunctionTreeVector<D> &inp,
                                  int max_iter) const {
    TreeBuilder<D> builder;
    CopyAdaptor<D> adaptor(inp, this->maxScale, 0);
    DefaultCalculator<D> calculator;
    builder.build(out, calculator, adaptor, max_iter);
    println(10, std::endl);
}

template class GridGenerator<1>;
template class GridGenerator<2>;
template class GridGenerator<3>;
