#include "MWOperator.h"
#include "FunctionTree.h"
#include "WaveletAdaptor.h"
#include "GridGenerator.h"
#include "OperApplicationCalculator.h"
#include "MultiResolutionAnalysis.h"
#include "Timer.h"

template<int D>
MWOperator<D>::MWOperator(const MultiResolutionAnalysis<D> &mra, double prec, int iter)
    : TreeBuilder<D>(mra, iter),
      apply_prec(prec) {
}

template<int D>
MWOperator<D>::~MWOperator() {
    this->clearAdaptor();
    this->clearOperator();
}

template<int D>
void MWOperator<D>::clearOperator() {
    for (int i = 0; i < this->oper.size(); i++) {
        if (this->oper[i] != 0) delete this->oper[i];
        this->oper.clear();
    }
}

template<int D>
FunctionTree<D>* MWOperator<D>::operator()(FunctionTree<D> &inp) {
    FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
    initializeGrid(*out, inp);
    (*this)(*out, inp);
    return out;
}

template<int D>
void MWOperator<D>::operator()(FunctionTree<D> &out, FunctionTree<D> &inp) {
    this->oper.calcBandWidths(this->apply_prec);
    this->adaptor = new WaveletAdaptor<D>(this->apply_prec);
    this->calculator = new OperApplicationCalculator<D>(this->oper, inp);
    this->build(out);
    this->clearCalculator();
    this->clearAdaptor();
    this->oper.clearBandWidths();

    Timer trans_t;
    trans_t.restart();
    out.mwTransform(TopDown, false); // add coarse scale contributions
    out.mwTransform(BottomUp);
    trans_t.stop();

    println(10, "Time transform      " << trans_t);
    println(10, std::endl);
}

/** Build grid based on analytic input function */
template<int D>
void MWOperator<D>::initializeGrid(FunctionTree<D> &out, FunctionTree<D> &inp) {
    Timer init_t;
    init_t.restart();
    GridGenerator<D> G(this->MRA);
    G(out, inp);
    init_t.stop();
    println(10, "Time initializing   " << init_t);
}

template<int D>
MultiResolutionAnalysis<2>* MWOperator<D>::getOperatorMRA() {
    const BoundingBox<D> &box = this->MRA.getWorldBox();
    const ScalingBasis &basis = this->MRA.getScalingBasis();

    int maxn = 0;
    for (int i = 0; i < D; i++) {
        if (box.size(i) > maxn) {
            maxn = box.size(i);
        }
    }
    int nbox[2] = { maxn, maxn};
    NodeIndex<2> idx(box.getScale());
    BoundingBox<2> oper_box(idx, nbox);
    return new MultiResolutionAnalysis<2>(oper_box, basis);
}

template class MWOperator<1>;
template class MWOperator<2>;
template class MWOperator<3>;
