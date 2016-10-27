#include "DerivativeOperator.h"
#include "DerivativeGenerator.h"

template<int D>
DerivativeOperator<D>::DerivativeOperator(int dir,
                                          const MultiResolutionAnalysis<D> &mra,
                                          double a,
                                          double b)
        : MWOperator<D>(mra, MachineZero),
          A(a),
          B(b) {
    if (this->A > MachineZero) NEEDS_TESTING;
    if (this->B > MachineZero) NEEDS_TESTING;
    initializeOperator();
    this->setApplyDir(dir);
}

template<int D>
DerivativeOperator<D>::~DerivativeOperator() {
    this->clearOperator();
}

template<int D>
void DerivativeOperator<D>::initializeOperator() {
    MultiResolutionAnalysis<2> *oper_mra = this->MRA.getOperatorMRA();
    DerivativeGenerator DG(this->MRA.getScalingBasis());

    OperatorTree *oper_comp = new OperatorTree(*oper_mra, MachineZero, MaxScale);
    DG(*oper_comp, this->A, this->B);
    this->oper.push_back(oper_comp);

    delete oper_mra;
}

template<int D>
void DerivativeOperator<D>::grad(FunctionTreeVector<D> &out,
                                 FunctionTree<D> &inp) {
    NOT_IMPLEMENTED_ABORT;
    /*
    GridGenerator<D> G(this->MRA);
    FunctionTreeVector<D> out;
    for (int d = 0; d < D; d++) {
        this->setApplyDir(d);
        FunctionTree<D> *out_d = G(inp);
        (*this)(*out_d, inp, 0);
        out.push_back(out_d);
    }
    return out;
    */
}

template<int D>
void DerivativeOperator<D>::div(FunctionTree<D> &out,
                                FunctionTreeVector<D> &inp) {
    NOT_IMPLEMENTED_ABORT;
    /*
    if (inp.size() != D) MSG_ERROR("Invalid dimension");

    MWAdder<D> add;
    GridGenerator<D> G;

    FunctionTreeVector<D> vec;
    for (int d = 0; d < D; d++) {
        this->setApplyDir(d);
        FunctionTree<D> *tmp_d = G(inp);
        (*this)(*tmp_d, *inp[d], 0);
        vec.push_back(tmp_d);
    }

    FunctionTree<D> *out = G(inp);
    add(*out, vec, 0);
    vec.clear(true);

    return out;
    */
}

template class DerivativeOperator<1>;
template class DerivativeOperator<2>;
template class DerivativeOperator<3>;
