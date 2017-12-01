#include "MWProjector.h"
#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "ProjectionCalculator.h"
#include "AnalyticFunction.h"
#include "FunctionTree.h"
#include "Timer.h"
#include "Printer.h"

template<int D>
void MWProjector<D>::operator()(FunctionTree<D> &out,
                                std::function<double (const double *r)> func,
                                int maxIter) const {
    AnalyticFunction<D> inp(func);
    (*this)(out, inp, maxIter);
}

template<int D>
void MWProjector<D>::operator()(FunctionTree<D> &out,
                                RepresentableFunction<D> &inp,
                                int maxIter) const {
    TreeBuilder<D> builder;
    WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
    ProjectionCalculator<D> calculator(inp);

    builder.build(out, calculator, adaptor, maxIter);

    Timer trans_t;
    out.mwTransform(BottomUp);
    out.calcSquareNorm();
    trans_t.stop();

    println(10, "Time transform      " << trans_t);
    println(10, std::endl);
}

template class MWProjector<1>;
template class MWProjector<2>;
template class MWProjector<3>;
